#ifndef _MGARD_RECOMPOSE_HPP
#define _MGARD_RECOMPOSE_HPP

#include <vector>
#include <cstdlib>

namespace MGARD{

using namespace std;

template <class T>
void recover_from_interpolant_difference(T const * data, T * coeff, size_t n, size_t stride, bool even){
	T const * prev = data;			// adjacent data in N_l (left)
	T const * next = data + stride;	// adjacent data in N_l (right)
	for(int i=0; i<n-1; i++){
		*coeff += (*prev + *next) / 2; 
		prev = next;
		next += stride;
		coeff += stride;
	}
	if(even){
		// n in N_l is even, deal with the last coefficient
		*coeff += *prev;
	}	
}

template <class T>
void subtract_correction(T * data, T const * correction, size_t n, size_t stride){
	T * data_pos = data;
	for(int i=0; i<n; i++){
		*data_pos -= correction[i];
		data_pos += stride;
	}
}

// recompose a level with n element and the given stride
// from a level with n/2 element
template <class T>
void recompose_level_1D(T * data, const vector<size_t>& levels, const vector<size_t>& strides, int current_level){
	size_t stride = strides[current_level];
	size_t n = levels[current_level];
	size_t next_stride = strides[current_level + 1];
	auto load_v = compute_load_vector(data + next_stride, n, stride);
	auto correction = compute_correction(load_v.data(), n, (T)stride);
	subtract_correction(data, correction.data(), n, stride);
	recover_from_interpolant_difference(data, data + next_stride, n, stride, !(n&1));
}

template <class T>
void recompose_level(T * data, const vector<size_t>& levels, const vector<size_t>& strides, int current_level){
	// decompose along each dimension
	recompose_level_1D(data, levels, strides, current_level);
}

template <class T>
void recompose(T * data, size_t n, size_t target_level){
	size_t stride = 1;
	vector<size_t> levels(target_level + 1);
	vector<size_t> strides(target_level + 1);
	for(int i=0; i<=target_level; i++){
		levels[target_level - i] = n;
		strides[target_level - i] = stride; 
		stride <<= 1;
		n = (n+1) >> 1;
	}
	for(int i=0; i<target_level; i++){
		recompose_level(data, levels, strides, i);
	}
}

}

#endif