#ifndef _MGARD_RECOMPOSE_HPP
#define _MGARD_RECOMPOSE_HPP

#include <vector>
#include <cstdlib>

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
	T * nodal_buffer = NULL;	// starting position for nodal values
	T * coeff_buffer = NULL;	// starting position for coefficients
	T * load_v_buffer = NULL;
	T * correction_buffer = NULL;

	// reorder the data to original order (insert coeffcients between nodal values)
	void data_reorder(T * data_pos, size_t n){
		const T * nodal_pos = nodal_buffer;
		const T * coeff_pos = coeff_buffer;
		T * cur_data_pos = data_pos;
		for(int i=0; i<n_coeff; i++){
			*(cur_data_pos++) = *(nodal_pos++);
			*(cur_data_pos++) = *(coeff_pos++);
		}
		*(cur_data_pos++) = *(nodal_pos++);
		if(!(n & 1)){
			// if even, the last coefficient equals to the interpolant
			// of the last two nodal values
			*cur_data_pos = (nodal_pos[-1] + nodal_pos[0]) / 2;
		}
	}

	void recover_from_interpolant_difference(){
		for(int i=0; i<n_coeff; i++){
			coeff_buffer[i] += (nodal_buffer[i] + nodal_buffer[i+1]) / 2; 
		}
	}

	void subtract_correction(){
		for(int i=0; i<n_nodal; i++){
			nodal_buffer[i] -= correction_buffer[i];
		}
	}

	// compute entries for load vector
	// for uniform decomposition only
	// @param h: stride of nodals in N_(l+1)
	// output in load_v_buffer
	void compute_load_vector(T h){
		T const * coeff = coeff_buffer;
		T ah = h * 0.25; // derived constant in the formula
		// first nodal value
		load_v_buffer[0] = coeff[0] * ah;
		// iterate through nodal values
		for(int i=1; i<n_coeff; i++){
			load_v_buffer[i] = (coeff[i-1] + coeff[i]) * ah;
		}
		// last nodal value
		load_v_buffer[n_coeff] = coeff[n_coeff-1] * ah;
		// if next n is even, load_v_buffer[n_nodal - 1] = 0
		if(n_nodal > n_coeff + 1) load_v_buffer[n_coeff + 1] = 0;
	}
	// compute correction on nodal value 
	// using Thomas algorithm for tridiagonal inverse
	// @param h: interval length
	void compute_correction(T h){
		size_t n = n_nodal;
		// Thomas algorithm for solving M_l x = load_v
		// forward pass
		// simplified algorithm
		T * d = load_v_buffer;
		vector<T> b(n, h*2/3);
		b[0] = h/3;
		b[n-1] = h/3;
		T c = h/6;
		for(int i=1; i<n; i++){
			auto w = c / b[i-1];
			b[i] = b[i] - w * c;
			d[i] = d[i] - w * d[i-1];
		}
		// backward pass
		correction_buffer[n-1] = d[n-1] / b[n-1];
		for(int i=n-2; i>=0; i--){
			correction_buffer[i] = (d[i] - c * correction_buffer[i+1]) / b[i];
		}
	}
	// recompose a level with n element and the given stride
	// from a level with n/2 element
	void recompose_level_1D(T * data_pos, size_t n, size_t stride, T h){
		cerr << n << endl;
		memcpy(data_buffer, data_pos, n*sizeof(T));
		n_nodal = (n >> 1) + 1;
		n_coeff = n - n_nodal;
		nodal_buffer = data_buffer;
		coeff_buffer = data_buffer + n_nodal;
		compute_load_vector(h);
		compute_correction(h);
		subtract_correction();
		recover_from_interpolant_difference();
		data_reorder(data_pos, n);
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