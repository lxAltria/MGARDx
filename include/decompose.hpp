#include <vector>
#include <cstdlib>

namespace MGARD{

using namespace std;

// compute the difference between original value 
// and interpolant (I - PI_l)Q_l
// overwrite the data in N_l \ N_(l-1) in place
template <class T>
void compute_interpolant_difference(T const * data, T * coeff, size_t n, size_t stride, bool even){
	T const * prev = data;					// adjacent data in N_l (left)
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
	// b[0] = h/3, b[1:n] = h/4	
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
	auto load_v = compute_load_vector(data + stride, next_n, next_stride);
	auto correction = compute_correction(load_v.data(), next_n, (T)next_stride);
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
