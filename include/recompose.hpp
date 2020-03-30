#ifndef _MGARD_RECOMPOSE_HPP
#define _MGARD_RECOMPOSE_HPP

#include <vector>
#include <cstdlib>
#include <cstring>
#include "utils.hpp"

namespace MGARD{

using namespace std;

template <class T>
class Recomposer{
public:
	Recomposer(){};
	~Recomposer(){
		if(data_buffer) free(data_buffer);
		if(correction_buffer) free(correction_buffer);	
		if(load_v_buffer) free(load_v_buffer);
	};
	void recompose(T * data_, const vector<size_t>& dims, size_t target_level){
		data = data_;
		vector<vector<size_t>> level_dims;
		// compute n_nodal in each level
		for(int i=0; i<dims.size(); i++){
			level_dims.push_back(vector<size_t>(target_level + 1));
			int n = dims[i];
			for(int j=0; j<=target_level; j++){
				level_dims[i][target_level - j] = n;
				n = (n >> 1) + 1;
			}
			cerr << "Ns: ";
			for(int j=0; j<=target_level; j++){
				cerr << level_dims[i][j] << " ";
			}
			cerr << endl;
		}
		size_t num_elements = 1;
		for(const auto& d:dims){
			num_elements *= d;
		}
		data_buffer_size = num_elements * sizeof(T);
		init(dims);
		size_t h = 1 << (target_level - 1);
		if(dims.size() == 1){
			for(int i=0; i<target_level; i++){
				recompose_level_1D(data, level_dims[0][i+1], h);
				h >>= 1;
			}
		}
		else if(dims.size() == 2){
			for(int i=0; i<target_level; i++){
				size_t n1 = level_dims[0][i+1];
				size_t n2 = level_dims[1][i+1];
				recompose_level_2D(data, n1, n2, (T)h, dims[1]);
				h >>= 1;
			}
		}
	}
private:
	unsigned int default_batch_size = 32;
	size_t data_buffer_size = 0;
	T * data = NULL;			// pointer to the original data
	T * data_buffer = NULL;		// buffer for reordered data
	T * load_v_buffer = NULL;
	T * correction_buffer = NULL;

	void init(const vector<size_t>& dims){
		size_t buffer_size = default_batch_size * (*max_element(dims.begin(), dims.end())) * sizeof(T);
		cerr << "buffer_size = " << buffer_size << endl;
		cerr << "data_buffer_size = " << data_buffer_size << endl;
		if(data_buffer) free(data_buffer);
		if(correction_buffer) free(correction_buffer);
		if(load_v_buffer) free(load_v_buffer);
		data_buffer = (T *) malloc(data_buffer_size);
		correction_buffer = (T *) malloc(buffer_size);
		load_v_buffer = (T *)malloc(buffer_size);
	}
	// reorder the data to original order (insert coeffcients between nodal values)
	void data_reverse_reorder_1D(T * data_pos, int n_nodal, int n_coeff, const T * nodal_buffer, const T * coeff_buffer){
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
	// reorder the data to recover the data order in next level
	/*
		oooxx		oooxx		oxoxo
		oooxx	(1)	xxxxx	(2)	xxxxo
		oooxx	=>	oooxx	=>	oxoxo
		xxxxx		xxxxx		xxxxx
		xxxxx		oooxx		oxoxo
	*/
	void data_reverse_reorder_2D(T * data_pos, size_t n1_nodal, size_t n1_coeff, size_t n2_nodal, size_t n2_coeff, size_t stride){
		size_t n1 = n1_nodal + n1_coeff;
		size_t n2 = n2_nodal + n2_coeff;
		T * cur_data_pos = data_pos;
		T * nodal_pos = data_buffer;
		T * coeff_pos = data_buffer + n2_nodal;
		// do reorder (1)
		// TODO: change to online processing for memory saving
		switch_rows_2D_by_buffer_reverse(data_pos, data_buffer, n1_nodal + n1_coeff, n2_nodal + n2_coeff, stride);
		// do reorder (2)
		for(int i=0; i<n1; i++){
			memcpy(data_buffer, cur_data_pos, n2 * sizeof(T));
			data_reverse_reorder_1D(cur_data_pos, n2_nodal, n2_coeff, nodal_pos, coeff_pos);
			cur_data_pos += stride;
		}
		if(!(n1 & 1)){
			// n1 is even, recover the coefficients
			cur_data_pos -= stride;
			for(int j=0; j<n2; j++){
				cur_data_pos[j] = (cur_data_pos[j] + cur_data_pos[-stride + j]) / 2;
			}
		}
	}
	void recover_from_interpolant_difference_1D(size_t n_coeff, const T * nodal_buffer, T * coeff_buffer){
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
	void recompose_level_1D(T * data_pos, size_t n, T h){
		cerr << n << endl;
		size_t n_nodal = (n >> 1) + 1;
		size_t n_coeff = n - n_nodal;
		memcpy(data_buffer, data_pos, n*sizeof(T));
		T * nodal_buffer = data_buffer;
		T * coeff_buffer = data_buffer + n_nodal;
		compute_load_vector(load_v_buffer, n_nodal, n_coeff, h, coeff_buffer);
		compute_correction(correction_buffer, n_nodal, h, load_v_buffer);
		subtract_correction(n_nodal, nodal_buffer);
		recover_from_interpolant_difference_1D(n_coeff, nodal_buffer, coeff_buffer);
		data_reverse_reorder_1D(data_pos, n_nodal, n_coeff, nodal_buffer, coeff_buffer);
	}
	// compute and subtract the corrections
	void compute_and_subtract_correction_2D(T * data_pos, size_t n1_nodal, size_t n1_coeff, size_t n2_nodal, size_t n2_coeff, T h, size_t stride){
		size_t n1 = n1_nodal + n1_coeff;
		size_t n2 = n2_nodal + n2_coeff;
		// compute horizontal correction
		T * nodal_pos = data_pos;
		const T * coeff_pos = data_pos + n2_nodal;
		// store horizontal corrections in the data_buffer
		T * correction_pos = data_buffer;
		for(int i=0; i<n1; i++){
			compute_load_vector(load_v_buffer, n2_nodal, n2_coeff, h, coeff_pos);
			compute_correction(correction_pos, n2_nodal, h, load_v_buffer);
			// subtract_correction(n2_nodal, nodal_pos);
			nodal_pos += stride, coeff_pos += stride;
			correction_pos += n2_nodal;
		}
		// compute vertical correction
		compute_and_apply_correction_2D_vertical(data_pos, n1, n2, h, stride, data_buffer, load_v_buffer, correction_buffer, default_batch_size, false);
	}
	// compute the difference between original value 
	// and interpolant (I - PI_l)Q_l for the coefficient rows in 2D
	// overwrite the data in N_l \ N_(l-1) in place
	// Note: interpolant difference in the nodal rows have already been computed
	void recover_from_interpolant_difference_2D_vertical(T * data_pos, size_t n1, size_t n2, size_t stride){
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
                *(nodal_coeff_pos++) += (nodal_pos[j] + nodal_pos[stride + j]) / 2;
                // coefficients in centers
                *(coeff_coeff_pos++) += (nodal_pos[j] + nodal_pos[j + 1] + nodal_pos[stride + j] + nodal_pos[stride + j + 1]) / 4;
            }
            // compute the last (or second last if n2 is even) nodal column
            *(nodal_coeff_pos ++) += (nodal_pos[n2_coeff] + nodal_pos[stride + n2_coeff]) / 2;
            if(even_n2){
                // compute the last nodal column
                *(nodal_coeff_pos ++) += (nodal_pos[n2_coeff + 1] + nodal_pos[stride + n2_coeff + 1]) / 2;
            }
		}
	}
	void recover_from_interpolant_difference_2D(T * data_pos, size_t n1_nodal, size_t n1_coeff, size_t n2_nodal, size_t n2_coeff, size_t stride){
		size_t n1 = n1_nodal + n1_coeff;
		size_t n2 = n2_nodal + n2_coeff;
		// compute horizontal difference
		const T * nodal_pos = data_pos;
		T * coeff_pos = data_pos + n2_nodal;
		for(int i=0; i<n1_nodal; i++){
			recover_from_interpolant_difference_1D(n2_coeff, nodal_pos, coeff_pos);
			nodal_pos += stride, coeff_pos += stride;
		}
		// compute vertical difference
		recover_from_interpolant_difference_2D_vertical(data_pos, n1, n2, stride);
	}	
	// recompose n1/2 x n2/2 data into fine level (n1 x n2)
	void recompose_level_2D(T * data_pos, size_t n1, size_t n2, T h, size_t stride){
		cerr << "recompose, h = " << h << endl; 
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		size_t n2_nodal = (n2 >> 1) + 1;
		size_t n2_coeff = n2 - n2_nodal;
		// if(h > 1) print(data, n1, n2, "data_after_correction");
		compute_and_subtract_correction_2D(data_pos, n1_nodal, n1_coeff, n2_nodal, n2_coeff, h, stride);
		// if(h > 1) print(data, n1, n2, "data_before_correction");
		recover_from_interpolant_difference_2D(data_pos, n1_nodal, n1_coeff, n2_nodal, n2_coeff, stride);
		data_reverse_reorder_2D(data_pos, n1_nodal, n1_coeff, n2_nodal, n2_coeff, stride);
		// if(h > 1) print(data, n1, n2, "data_entry");
	}
};

}

#endif