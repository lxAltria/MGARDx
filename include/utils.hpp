#ifndef _MGARD_UTILS_HPP
#define _MGARD_UTILS_HPP

#include <vector>
#include <cstdlib>

namespace MGARD{

using namespace std;

// compute entries for load vector
// for uniform decomposition only
// @param n_nodal: number of nodal points (N_l)
// @param n_coeff: number of coefficients (N_(l+1) - N_l)
// @param coeff: starting location of N_(l+1) \ N_l
// @param h: stride of nodals in N_(l+1)
// output in load_v_buffer
template <class T>
void compute_load_vector(size_t n_nodal, size_t n_coeff, T const * coeff, T h, T * load_v_buffer){
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
// @param data: starting location
// @param load_v: load vector
// @param h: interval length
// @param n: number of points in current level
// output in result
template <class T>
void compute_correction(T * load_v_buffer, size_t n, T h, T * result){
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
	result[n-1] = d[n-1] / b[n-1];
	for(int i=n-2; i>=0; i--){
		result[i] = (d[i] - c * result[i+1]) / b[i];
	}
}

}

#endif