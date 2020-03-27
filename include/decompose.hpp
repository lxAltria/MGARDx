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
		data_buffer_size = stride * sizeof(T);
		init();
		if(dims.size() == 1){
			for(int i=0; i<target_level; i++){
				decompose_level();
				for(int i=0; i<dims.size(); i++){
					dims[i] = (dims[i] >> 1) + 1;
					h_l[i] <<= 1;
				}
			}
		}
		else if(dims.size() == 2){
			size_t h = 1;
			size_t n1 = dims[0];
			size_t n2 = dims[1];
			for(int i=0; i<target_level; i++){
				decompose_level_2D(data, n1, n2, (T)h, n2);
				n1 = (n1 >> 1) + 1;
				n2 = (n2 >> 1) + 1;
				h <<= 1;
			}
		}
	}

private:
	unsigned int default_batch_size = 1;
	size_t data_buffer_size = 0;
	vector<size_t> dims;		// dims for each level MUTABLE
	vector<size_t> h_l;			// interval length for each dimension MUTABLE
	vector<size_t> strides;		// strides for each dimension CONSTANT
	T * data = NULL;			// pointer to the original data
	T * data_buffer = NULL;		// buffer for reordered data
	T * load_v_buffer = NULL;
	T * correction_buffer = NULL;

	void init(){
		size_t buffer_size = default_batch_size * (*max_element(dims.begin(), dims.end())) * sizeof(T);
		cerr << "buffer_size = " << buffer_size << endl;
		if(data_buffer) free(data_buffer);
		if(correction_buffer) free(correction_buffer);
		if(load_v_buffer) free(load_v_buffer);
		data_buffer = (T *) malloc(data_buffer_size);
		correction_buffer = (T *) malloc(buffer_size);
		load_v_buffer = (T *)malloc(buffer_size);
	}
	// reorder the data to put all the coefficient to the back
	void data_reorder_1D(int n_nodal, int n_coeff, T const * data_pos, T * nodal_buffer, T * coeff_buffer){
		T * nodal_pos = nodal_buffer;
		T * coeff_pos = coeff_buffer;
		T const * cur_data_pos = data_pos;
		for(int i=0; i<n_coeff; i++){
			*(nodal_pos++) = *(cur_data_pos++);
			*(coeff_pos++) = *(cur_data_pos++);
		}
		*(nodal_pos++) = *(cur_data_pos++);
		if(n_nodal == n_coeff + 2){
			// if even, add a nodal value such that the interpolant
			// of the last two nodal values equal to the last coefficient
			*nodal_pos = 2*cur_data_pos[0] - nodal_pos[-1];
		}
	}
	// compute the difference between original value 
	// and interpolant (I - PI_l)Q_l
	// overwrite the data in N_l \ N_(l-1) in place
	void compute_interpolant_difference(size_t n_coeff, const T * nodal_buffer, T * coeff_buffer){
		for(int i=0; i<n_coeff; i++){
			coeff_buffer[i] -= (nodal_buffer[i] + nodal_buffer[i+1]) / 2; 
		}
	}
	void add_correction(size_t n_nodal, T * nodal_buffer){
		for(int i=0; i<n_nodal; i++){
			nodal_buffer[i] += correction_buffer[i];
		}
	}
	// decompose a level with n element and the given stride
	// to a level with n/2 element
	void decompose_level_1D(T * data_pos, size_t n, T h){
		size_t n_nodal = (n >> 1) + 1;
		size_t n_coeff = n - n_nodal;
		T * nodal_buffer = data_buffer;
		T * coeff_buffer = data_buffer + n_nodal;
		data_reorder_1D(n_nodal, n_coeff, data_pos, nodal_buffer, coeff_buffer);
		compute_interpolant_difference(n_coeff, nodal_buffer, coeff_buffer);
		compute_load_vector(n_nodal, n_coeff, h, coeff_buffer, load_v_buffer);
		compute_correction(n_nodal, h, load_v_buffer, correction_buffer);
		add_correction(n_nodal, nodal_buffer);
		memcpy(data_pos, data_buffer, n*sizeof(T));
	}
	// switch the rows in the data for coherent memory access
	// o: nodal data, x: coefficient data
	/*
		oooxx		oooxx
		xxxxx		oooxx
		oooxx	=>	oooxx
		xxxxx		xxxxx
		oooxx		xxxxx
	*/
	void switch_rows_2D_by_buffer(T * data_pos, size_t n1, size_t n2, size_t stride){
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		T * nodal_data_buffer = data_buffer;
		T * coeff_data_buffer = data_buffer + n1_nodal * n2;
		T * cur_data_pos = data_pos + stride;
		for(int i=1; i<n1_coeff; i++){
			// copy coefficient rows
			memcpy(coeff_data_buffer + i * n2, cur_data_pos, n2 * sizeof(T));
			cur_data_pos += stride;
			// copy nodal rows
			memcpy(nodal_data_buffer + i * n2, cur_data_pos, n2 * sizeof(T));
			cur_data_pos += stride;
		}
		if(!(n1&1)){
			// n1 is even, move the last nodal row
			memcpy(coeff_data_buffer - stride, cur_data_pos, n2 * sizeof(T));
		}
		// copy data back
		cur_data_pos = data_pos + stride;
		for(int i=1; i<n1; i++){
			memcpy(cur_data_pos, data_buffer + i * n2, n2 * sizeof(T));
			cur_data_pos += stride;
		}
	}
	// compute the difference between original value 
	// and interpolant (I - PI_l)Q_l for the coefficient rows in 2D
	// overwrite the data in N_l \ N_(l-1) in place
	// Note: interpolant difference in the nodal rows have already been computed
	void compute_interpolant_difference_2D_vertical(T * data_pos, size_t n1, size_t n2, size_t stride){
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		size_t n2_nodal = (n2 >> 1) + 1;
		size_t n2_coeff = n2 - n2_nodal;
		bool even_n2 = !(n2 & 1);
		T * n1_nodal_data = data_pos;
		T * n1_coeff_data = data_pos + n1_nodal * stride;
		for(int i=0; i<n1_coeff; i++){
            const T * nodal_pos = n1_nodal_data + i * stride;
            T * coeff_pos = n1_coeff_data + i * stride;
            // TODO: optimize average computation
            T * nodal_coeff_pos = coeff_pos;	// coeffcients in nodal rows
            T * coeff_coeff_pos = coeff_pos + n2_nodal;	// coefficients in coeffcients rows
            for(int j=0; j<n2_coeff; j++){
                // coefficients in nodal columns
                *(nodal_coeff_pos++) -= (nodal_pos[j] + nodal_pos[stride + j]) / 2;
                // coefficients in centers
                *(coeff_coeff_pos++) -= (nodal_pos[j] + nodal_pos[j + 1] + nodal_pos[stride + j] + nodal_pos[stride + j + 1]) / 4;
            }
            // compute the last (or second last if n2 is even) nodal column
            *(nodal_coeff_pos ++) -= (nodal_pos[n2_coeff] + nodal_pos[stride + n2_coeff]) / 2;
            if(even_n2){
                // compute the last nodal column
                *(nodal_coeff_pos ++) -= (nodal_pos[n2_coeff + 1] + nodal_pos[stride + n2_coeff + 1]) / 2;
            }
		}
	}
	// compute entries for load vector in vertical (non-contiguous) direction
	// for uniform decomposition only
	// @param n: dimensions
	// @param stride: stride for vertical adjacent data
	// @param h: stride of nodals in N_(l+1)
	// @param batchsize: number of columns to be computed together
	// output in load_v_buffer
	void compute_load_vector_vertical(const T * coeff_buffer, size_t n1_nodal, size_t n1_coeff, size_t stride, T h, int batchsize){
		T ah = h * 0.25; // derived constant in the formula
		T const * coeff_pos = coeff_buffer;
		T * load_v_pos = load_v_buffer;
		// first nodal value
		for(int j=0; j<batchsize; j++){
			load_v_pos[j] = coeff_pos[j] * ah;
		}
		load_v_pos += batchsize;
		coeff_pos += stride;
		for(int i=1; i<n1_coeff; i++){
			for(int j=0; j<batchsize; j++){
				load_v_pos[j] = (coeff_pos[-stride + j] + coeff_pos[j]) * ah;
			}
			load_v_pos += batchsize;
			coeff_pos += stride;
		}
		// last nodal value
		for(int j=0; j<batchsize; j++){
			load_v_pos[j] = coeff_pos[-stride + j] * ah;
		}
		// if next n is even, load_v_buffer[n_nodal - 1] = 0
		if(n1_nodal == n1_coeff + 2){
			for(int j=0; j<batchsize; j++){
				load_v_pos[j] = 0;
			}
		}
	}
	// compute correction on nodal value 
	// using Thomas algorithm for tridiagonal inverse
	// @param h: interval length
	void compute_correction_batched(T h, const T * b, const T * w, size_t n_nodal, int batchsize){
		size_t n = n_nodal;
		T c = h/6;
		// Thomas algorithm for solving M_l x = load_v
		// forward pass
		// simplified algorithm
		// b[:], w[:] are precomputed
		T * load_v_pos = load_v_buffer + batchsize;
		for(int i=1; i<n; i++){
			for(int j=0; j<batchsize; j++){
				load_v_pos[j] -= - w[i] * load_v_pos[-batchsize + j];
			}
			load_v_pos += batchsize;
		}
		// backward pass
		T * correction_pos = correction_buffer + (n - 1) * batchsize;
		for(int j=0; j<batchsize; j++){
			correction_pos[j] = load_v_pos[j] / b[n - 1];
		}
		correction_pos -= batchsize;
		load_v_pos -= batchsize;
		for(int i=n-2; i>=0; i--){
			for(int j=0; j<batchsize; j++){
				correction_pos[j] = (load_v_pos[j] - c * correction_pos[batchsize + j]) / b[i];
			}
			correction_pos -= batchsize;
			load_v_pos -= batchsize;
		}
	}
	void add_correction_batched(T * nodal_pos, int n_nodal, int stride, int batchsize){
		const T * correction_pos = correction_buffer;
		for(int i=0; i<n_nodal; i++){
			for(int j=0; j<batchsize; j++){
				nodal_pos[j] += correction_pos[j];
			}
			nodal_pos += stride;
			correction_pos += batchsize;
		}
	}
	// decompose 2D in the vertical (non-contiguous) dimension
	void decompose_level_2D_vertical(T * data_pos, size_t n1, size_t n2, T h, size_t stride){
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		size_t n2_nodal = (n2 >> 1) + 1;
		size_t n2_coeff = n2 - n2_nodal;
		compute_interpolant_difference_2D_vertical(data_pos, n1, n2, stride);
		vector<T> b(n1_nodal, h*2/3);
		vector<T> w(n1_nodal, 0);
		b[0] = h/3,	b[n1_nodal - 1] = h/3;
		T c = h/6;
		for(int i=1; i<n1_nodal; i++){
			w[i] = c / b[i-1];
			b[i] = b[i] - w[i] * c;
		}
		int batchsize = default_batch_size;
		int num_batches = (n2_coeff - 1) / batchsize;
		T * nodal_pos = data_pos;
		T * coeff_pos = data_pos + n1_nodal * stride;
		for(int i=0; i<num_batches; i++){
			compute_load_vector_vertical(coeff_pos, n1_nodal, n1_coeff, stride, h, batchsize);
			compute_correction_batched(h, b.data(), w.data(), n1_nodal, batchsize);
			add_correction_batched(nodal_pos, n1_nodal, stride, batchsize);
			nodal_pos += batchsize, coeff_pos += batchsize;
		}
		if(n2_coeff - batchsize * num_batches > 0){
			batchsize = n2_coeff - batchsize * num_batches;
			compute_load_vector_vertical(coeff_pos, n1_nodal, n1_coeff, stride, h, batchsize);
			compute_correction_batched(h, b.data(), w.data(), n1_nodal, batchsize);
			add_correction_batched(nodal_pos, n1_nodal, stride, batchsize);
		}
	}
	// horizontal sweep: decompose both coefficient and nodal rows
	void decompose_level_2D_horizontal(T * data_pos, size_t n1, size_t n2, T h, size_t stride){
		if(!(n1&1)){
			// n1 is even, change the last row into nodal row
			T * cur_data_pos = data_pos + (n1 - 1) * stride;
			for(int j=0; j<n2; j++){
				cur_data_pos[j] = 2 * cur_data_pos[j] - cur_data_pos[- stride + j];
			}
		}
		T * cur_data_pos = data_pos;
		for(int i=0; i<n1; i++){
			decompose_level_1D(cur_data_pos, n2, h);
			cur_data_pos += stride;
		}
	}
	// decompose n1 x n2 data into coarse level (n1/2 x n2/2)
	void decompose_level_2D(T * data_pos, size_t n1, size_t n2, T h, size_t stride){
		// decompose along dimension n2
		/*
			oxoxo		oooxx
			xxxxx		xxxxx
			oxoxo	=>	oooxx
			xxxxx		xxxxx
			oxoxo		oooxx
		*/
		decompose_level_2D_horizontal(data_pos, n1, n2, h, stride);
		// switch rows
		/*
			oooxx		oooxx
			xxxxx		oooxx
			oooxx	=>	oooxx
			xxxxx		xxxxx
			oooxx		xxxxx
		*/
		if(n1 * n2 <= data_buffer_size){
			switch_rows_2D_by_buffer(data_pos, n1, n2, stride);
		}
		else{
			// TODO: memcpy row by row
			cerr << "Not supported\n" << endl;
			exit(0);
		}
		// decompose along n2
		decompose_level_2D_vertical(data_pos, n1, n2, h, stride);
	}

	void decompose_level(){
		// decompose along each dimension
		for(int i=0; i<dims.size(); i++){
			decompose_level_1D(data, dims[i], (T) h_l[i]);
		}
	}
};


}

#endif