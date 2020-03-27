#ifndef _MGARD_RECOMPOSE_HPP
#define _MGARD_RECOMPOSE_HPP

#include <vector>
#include <cstdlib>
#include "utils.hpp"

namespace MGARD{

using namespace std;

template <class T>
class Recomposer{
public:
	Recomposer(){
		dims = vector<vector<size_t>>();
		h_l = vector<vector<size_t>>();
	};
	~Recomposer(){
		if(data_buffer) free(data_buffer);
		if(correction_buffer) free(correction_buffer);	
		if(load_v_buffer) free(load_v_buffer);
	};
	void recompose(T * data_, const vector<size_t>& dims_, size_t target_level){
		// compute n_nodal in each level
		for(int i=0; i<dims_.size(); i++){
			dims.push_back(vector<size_t>(target_level + 1));
			h_l.push_back(vector<size_t>(target_level + 1));
			int n = dims_[i];
			size_t h = 1;
			for(int j=0; j<=target_level; j++){
				dims[i][target_level - j] = n;
				h_l[i][target_level - j] = h; 
				n = (n >> 1) + 1;
				h <<= 1;
			}
			cerr << "Ns: ";
			for(int j=0; j<=target_level; j++){
				cerr << dims[i][j] << " ";
			}
			cerr << endl;
		}
		// compute constant strides of each dimension
		size_t stride = 1;
		strides.reserve(dims_.size());
		for(int i=0; i<dims_.size(); i++){
			strides[i] = stride;
			stride *= dims_[i];
		}
		// init buffers 
		size_t buffer_size = batchsize * (*max_element(dims_.begin(), dims_.end())) * sizeof(T);
		if(data_buffer) free(data_buffer);
		if(correction_buffer) free(correction_buffer);
		if(load_v_buffer) free(load_v_buffer);
		data_buffer = (T *) malloc(buffer_size);
		correction_buffer = (T *) malloc(buffer_size);
		load_v_buffer = (T *)malloc(buffer_size);
		// recompose
		data = data_;		
		for(int i=0; i<target_level; i++){
			cerr << "recompose level " << i << endl;
			recompose_level(i);
		}
	}
private:
	size_t n_nodal = 0;
	size_t n_coeff = 0;
	unsigned int batchsize = 1;
	vector<vector<size_t>> dims;// dims for each level
	vector<vector<size_t>> h_l;	// interval length for each dimension
	vector<size_t> strides;		// strides for each dimension CONSTANT
	T * data = NULL;			// pointer to the original data
	T * data_buffer = NULL;		// buffer for reordered data
	T * load_v_buffer = NULL;
	T * correction_buffer = NULL;

	// reorder the data to original order (insert coeffcients between nodal values)
	void data_reverse_reorder_1D(int n_nodal, int n_coeff, T * data_pos, const T * nodal_buffer, const T * coeff_buffer){
		const T * nodal_pos = nodal_buffer;
		const T * coeff_pos = coeff_buffer;
		T * cur_data_pos = data_pos;
		for(int i=0; i<n_coeff; i++){
			*(cur_data_pos++) = *(nodal_pos++);
			*(cur_data_pos++) = *(coeff_pos++);
		}
		*(cur_data_pos++) = *(nodal_pos++);
		if(n_nodal == n_coeff + 2){
			// if even, the last coefficient equals to the interpolant
			// of the last two nodal values
			*cur_data_pos = (nodal_pos[-1] + nodal_pos[0]) / 2;
		}
	}
	void recover_from_interpolant_difference(size_t n_coeff, const T * nodal_buffer, T * coeff_buffer){
		for(int i=0; i<n_coeff; i++){
			coeff_buffer[i] += (nodal_buffer[i] + nodal_buffer[i+1]) / 2; 
		}
	}
	void subtract_correction(size_t n_nodal, T * nodal_buffer){
		for(int i=0; i<n_nodal; i++){
			nodal_buffer[i] -= correction_buffer[i];
		}
	}
	// recompose a level with n element and the given stride
	// from a level with n/2 element
	void recompose_level_1D(T * data_pos, size_t n, size_t stride, T h){
		cerr << n << endl;
		n_nodal = (n >> 1) + 1;
		n_coeff = n - n_nodal;
		memcpy(data_buffer, data_pos, n*sizeof(T));
		T * nodal_buffer = data_buffer;
		T * coeff_buffer = data_buffer + n_nodal;
		compute_load_vector(n_nodal, n_coeff, h, coeff_buffer, load_v_buffer);
		compute_correction(n_nodal, h, load_v_buffer, correction_buffer);
		subtract_correction(n_nodal, nodal_buffer);
		recover_from_interpolant_difference(n_coeff, nodal_buffer, coeff_buffer);
		data_reverse_reorder_1D(n_nodal, n_coeff, data_pos, nodal_buffer, coeff_buffer);
	}
	void recompose_level(size_t level){
		// decompose along each dimension
		for(int i=0; i<dims.size(); i++){
			recompose_level_1D(data, dims[i][level+1], strides[i], (T) h_l[i][level+1]);
		}
	}	
};

}

#endif