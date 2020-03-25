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

	// decompose a level with n element and the given stride
	// to a level with n/2 element
	void decompose_level_1D(T * data_pos, size_t n, size_t stride, T h){
		n_nodal = (n >> 1) + 1;
		n_coeff = n - n_nodal;
		nodal_buffer = data_buffer;
		coeff_buffer = data_buffer + n_nodal;
		data_reorder(data_pos, n);
		compute_interpolant_difference();
		compute_load_vector(h);
		compute_correction(h);
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