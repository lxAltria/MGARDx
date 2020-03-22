#ifndef _MGARD_UTILS_HPP
#define _MGARD_UTILS_HPP

#include <vector>
#include <cstdlib>

namespace MGARD{

using namespace std;

// compute entries for load vector
// for uniform decomposition only
// @param coeff: starting location of N_(l+1) \ N_l
// @param n: number of points in current level
// @param stride: stride in current level (N_l)
template <class T>
vector<T> compute_load_vector(T const * coeff, size_t n, size_t stride){
	vector<T> load_v_buffer(n);
	T const * prev = coeff;
	T const * next = coeff + stride;
	T ah = stride * 0.125; // derived constant in the formula
	// first nodal value
	load_v_buffer[0] = *prev * ah;
	// iterate through nodal values
	for(int i=1; i<n-1; i++){
		load_v_buffer[i] = (*prev + *next) * ah;
		prev = next;
		next += stride;
	}
	// last nodal value
	load_v_buffer[n-1] = *prev * ah;
	return load_v_buffer;
}

// compute correction on nodal value 
// using Thomas algorithm for tridiagonal inverse
// @param data: starting location
// @param load_v: load vector
// @param h: interval length
// @param n: number of points in current level
// output in correction_buffer
template <class T>
vector<T> compute_correction(T * load_v_buffer, size_t n, T h){
	// Thomas algorithm for solving M_l x = load_v
	// forward pass
	// simplified algorithm
	T * d = load_v_buffer;
	vector<T> b(n, h/3);
	T c = h/6;
	for(int i=1; i<n; i++){
		auto w = c / b[i-1];
		b[i] = b[i] - w * c;
		d[i] = d[i] - w * d[i-1];
	}
	// backward pass
	vector<T> result(n);
	result[n-1] = d[n-1] / b[n-1];
	for(int i=n-2; i>=0; i--){
		result[i] = (d[i] - c * result[i+1]) / b[i];
	}
	return result;
}

}

#endif