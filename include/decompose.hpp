#ifndef _MGARD_DECOMPOSE_HPP
#define _MGARD_DECOMPOSE_HPP

#include <vector>
#include <cstdlib>
#include "utils.hpp"

namespace MGARD{

using namespace std;

// compute the difference between original value 
// and interpolant (I - PI_l)Q_l
// overwrite the data in N_l \ N_(l-1) in place
template <class T>
void compute_interpolant_difference(T const * data, T * coeff, size_t n, size_t stride, bool even){
	T const * prev = data;			// adjacent data in N_l (left)
	T const * next = data + stride;	// adjacent data in N_l (right)
	for(int i=0; i<n-1; i++){
		*coeff -= (*prev + *next) / 2; 
		prev = next;
		next += stride;
		coeff += stride;
	}
	if(even){
		// n in N_l is even, deal with the last coefficient
		*coeff -= *prev;
	}	
}

template <class T>
void add_correction(T * data, T const * correction, size_t n, size_t stride){
	T * data_pos = data;
	for(int i=0; i<n; i++){
		*data_pos += correction[i];
		data_pos += stride;
	}
}

// decompose a level with n element and the given stride
// to a level with n/2 element
template <class T>
void decompose_level_1D(T * data, size_t n, size_t stride){
	size_t next_n = (n+1) >> 1;
	size_t next_stride = stride << 1;
	compute_interpolant_difference(data, data + stride, next_n, next_stride, !(n&1));
	auto load_v = compute_load_vector(data + stride, next_n, next_stride, !(n&1));
	auto correction = compute_correction(load_v.data(), next_n, (T)next_stride, !(n&1));
	add_correction(data, correction.data(), next_n, next_stride);
}

template <class T>
void decompose_level(T * data, size_t n, size_t stride){
	// decompose along each dimension
	decompose_level_1D(data, n, stride);
}

template <class T>
void decompose(T * data, size_t n, size_t target_level){
	size_t stride = 1;
	for(int i=0; i<target_level; i++){
		decompose_level(data, n, stride);
		stride <<= 1;
		n = (n+1) >> 1;
	}
}

}

#endif