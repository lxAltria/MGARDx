#ifndef _MGARD_DECOMPOSE_HPP
#define _MGARD_DECOMPOSE_HPP

#include <vector>
#include <cstdlib>
#include <algorithm>
#include "utils.hpp"

namespace MGARD{

using namespace std;

template <class T>
class Decomposer{
public:
	Decomposer(){};
	~Decomposer(){
		if(data_buffer) free(data_buffer);
		if(correction_buffer) free(correction_buffer);	
		if(load_v_buffer) free(load_v_buffer);
	};
	void decompose(T * data_, const vector<size_t>& dims_, size_t target_level){
		data = data_;
		dims = dims_;
		h_l.reserve(dims.size());
		for(int i=0; i<dims.size(); i++){
			h_l[i] = 1;
		}
		strides.reserve(dims.size());
		size_t stride = 1;
		for(int i=0; i<dims_.size(); i++){
			strides[i] = stride;
			stride *= dims_[i];
		}
		init();
		for(int i=0; i<target_level; i++){
			decompose_level();
			for(int i=0; i<dims.size(); i++){
				dims[i] = (dims[i] >> 1) + 1;
				h_l[i] <<= 1;
			}
		}
	}

private:
	size_t n_nodal = 0;
	size_t n_coeff = 0;
	unsigned int batchsize = 1;
	vector<size_t> dims;		// dims for each level MUTABLE
	vector<size_t> h_l;			// interval length for each dimension MUTABLE
	vector<size_t> strides;		// strides for each dimension CONSTANT
	T * data = NULL;			// pointer to the original data
	T * data_buffer = NULL;		// buffer for reordered data
	T * nodal_buffer = NULL;	// starting position for nodal values
	T * coeff_buffer = NULL;	// starting position for coefficients
	T * load_v_buffer = NULL;
	T * correction_buffer = NULL;

	void init(){
		size_t buffer_size = batchsize * (*max_element(dims.begin(), dims.end())) * sizeof(T);
		if(data_buffer) free(data_buffer);
		if(correction_buffer) free(correction_buffer);
		if(load_v_buffer) free(load_v_buffer);
		data_buffer = (T *) malloc(buffer_size);
		correction_buffer = (T *) malloc(buffer_size);
		load_v_buffer = (T *)malloc(buffer_size);
	}

	// reorder the data to put all the coefficient to the back
	void data_reorder(T const * data_pos, size_t n){
		T * nodal_pos = nodal_buffer;
		T * coeff_pos = coeff_buffer;
		T const * cur_data_pos = data_pos;
		for(int i=0; i<n_coeff; i++){
			*(nodal_pos++) = *(cur_data_pos++);
			*(coeff_pos++) = *(cur_data_pos++);
		}
		*(nodal_pos++) = *(cur_data_pos++);
		if(!(n & 1)){
			// if even, add a nodal value such that the interpolant
			// of the last two nodal values equal to the last coefficient
			*nodal_pos = 2*cur_data_pos[0] - nodal_pos[-1];
		}
	}

	// compute the difference between original value 
	// and interpolant (I - PI_l)Q_l
	// overwrite the data in N_l \ N_(l-1) in place
	void compute_interpolant_difference(){
		for(int i=0; i<n_coeff; i++){
			coeff_buffer[i] -= (nodal_buffer[i] + nodal_buffer[i+1]) / 2; 
		}
	}

	void add_correction(){
		for(int i=0; i<n_nodal; i++){
			nodal_buffer[i] += correction_buffer[i];
		}
	}

	// decompose a level with n element and the given stride
	// to a level with n/2 element
	void decompose_level_1D(T * data_pos, size_t n, size_t stride, T h){
		n_nodal = (n >> 1) + 1;
		n_coeff = n - n_nodal;
		nodal_buffer = data_buffer;
		coeff_buffer = data_buffer + n_nodal;
		data_reorder(data_pos, n);
		compute_interpolant_difference();
		compute_load_vector(n_nodal, n_coeff, coeff_buffer, h, load_v_buffer);
		compute_correction(load_v_buffer, n_nodal, h, correction_buffer);
		add_correction();
		memcpy(data_pos, data_buffer, n*sizeof(T));
	}

	void decompose_level(){
		// decompose along each dimension
		for(int i=0; i<dims.size(); i++){
			decompose_level_1D(data, dims[i], strides[i], (T) h_l[i]);
		}
	}
};


}

#endif