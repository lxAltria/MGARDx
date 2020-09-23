#ifndef _MGARD_DECOMPOSE_HPP
#define _MGARD_DECOMPOSE_HPP

#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include "adaptive.hpp"
#include "utils.hpp"
#include "sz_compress_3d.hpp"

namespace MGARD{

using namespace std;

template <class T>
class Decomposer{
public:
	Decomposer(bool use_sz_=true){
            use_sz = use_sz_;
        };
	~Decomposer(){
		if(data_buffer) free(data_buffer);
		if(correction_buffer) free(correction_buffer);	
		if(load_v_buffer) free(load_v_buffer);
	};
	unsigned char * compress(T * data_, const vector<size_t>& dims, size_t target_level, double eb, size_t& compressed_size){
        error_bound = eb;
		target_level = decompose(data_, dims, target_level);
		size_t num_elements = 1;
		for(const auto& d:dims){
			num_elements *= d;
		}
        // writefile(string("decomposed.dat").c_str(), data, num_elements);
		auto result = quantize_and_encoding(dims, num_elements, dims.size(), eb, target_level, compressed_size);
		return result;
	}
    // return levels
	int decompose(T * data_, const vector<size_t>& dims, size_t target_level){
		data = data_;
		size_t num_elements = 1;
		for(const auto& d:dims){
			num_elements *= d;
		}
		data_buffer_size = num_elements * sizeof(T);
        int max_level = log2(*min_element(dims.begin(), dims.end()));
        if(target_level > max_level) target_level = max_level;
		init(dims);
        current_dims.resize(dims.size());
		if(dims.size() == 1){
			size_t h = 1;
			size_t n = dims[0];
			for(int i=0; i<target_level; i++){
				decompose_level_1D(data, n, h);
				n = (n >> 1) + 1;
				h <<= 1;
			}
		}
		else if(dims.size() == 2){
			size_t h = 1;
			size_t n1 = dims[0];
			size_t n2 = dims[1];
			for(int i=0; i<target_level; i++){
				decompose_level_2D(data, n1, n2, (T)h, dims[1]);
				n1 = (n1 >> 1) + 1;
				n2 = (n2 >> 1) + 1;
				h <<= 1;
			}
		}
		else if(dims.size() == 3){
			size_t h = 1;
			size_t n1 = dims[0];
			size_t n2 = dims[1];
			size_t n3 = dims[2];
            double C2 = 1 + 3*sqrt(3)/4;
            double c = sqrt(8);
			for(int i=0; i<target_level; i++){
                double cc = (1 - c) / (1 - pow(c, i + 1));
                double eb = cc * error_bound / C2;
                if(switch_to_lorenzo(data, n1, n2, n3, dims[1] * dims[2], dims[2], eb)){
                    cout << "switch to SZ (lorenzo) at level " << i << endl;
                    return i;
                }
				decompose_level_3D(data, n1, n2, n3, (T)h, dims[1] * dims[2], dims[2]);
				n1 = (n1 >> 1) + 1;
				n2 = (n2 >> 1) + 1;
				n3 = (n3 >> 1) + 1;
				h <<= 1;
                current_dims[0] = n1;
                current_dims[1] = n2;
                current_dims[2] = n3;
                // cerr << current_dims[0] << " " << current_dims[1] << " " << current_dims[2] << " " << endl;
			}
		}
        return target_level;
	}

private:
    double error_bound = 1e-6;
	unsigned int default_batch_size = 32;
	size_t data_buffer_size = 0;
    bool use_sz = true;
	T * data = NULL;			// pointer to the original data
	T * data_buffer = NULL;		// buffer for reordered data
	T * load_v_buffer = NULL;
	T * correction_buffer = NULL;
    vector<size_t> current_dims;

    int quantize_level(T * data, const vector<size_t>& dims, const vector<size_t>& coarse_dims, const vector<size_t>& fine_dims, double level_eb, int quant_radius, vector<int>& quant_inds, int offset, unsigned char *& compressed_data_pos){
        auto quantizer = SZ::LinearQuantizer<T>(level_eb, quant_radius);
        int start_offset = offset;
        for(int i=0; i<fine_dims[0]; i++){
            for(int j=0; j<fine_dims[1]; j++){
                for(int k=0; k<fine_dims[2]; k++){
                    if((i < coarse_dims[0]) && (j < coarse_dims[1]) && (k < coarse_dims[2])){
                        continue;
                    }
                    auto tmp = data[i * dims[1] * dims[2] + j * dims[2] + k];
                    quant_inds[offset ++] = quantizer.quantize_and_overwrite(tmp, 0);
                }
            }
        }
        // record quantizer
        quantizer.save(compressed_data_pos);
        return offset - start_offset;
    }

	unsigned char * quantize_and_encoding(const vector<size_t>& dims, size_t num_elements, int n_dims, double eb, int target_level, size_t& compressed_size){
        size_t num_nodal_elements = 1;
        for(const auto& d: current_dims){
            num_nodal_elements *= d;
        }
        unsigned char * compressed = (unsigned char *) malloc(num_elements * sizeof(T));
        unsigned char * compressed_data_pos = compressed;
        // record target level
        *reinterpret_cast<size_t*>(compressed_data_pos) = target_level;
        compressed_data_pos += sizeof(size_t);
        if(target_level == 0){
            // use sz directly
            // record eb
            size_t sz_compressed_size = 0;
            auto sz_compressed = sz_compress_3d(data, dims[0], dims[1], dims[2], eb/4, sz_compressed_size, 6);
            // record sz compressed data
            *reinterpret_cast<size_t*>(compressed_data_pos) = sz_compressed_size;
            compressed_data_pos += sizeof(size_t);
            cout << "sz compress position = " << compressed_data_pos - compressed << endl;
            cout << "sz compressed size = " << sz_compressed_size << endl;
            memcpy(compressed_data_pos, sz_compressed, sz_compressed_size);
            compressed_data_pos += sz_compressed_size;
            free(sz_compressed);            
        }
        else{
            *reinterpret_cast<unsigned char*>(compressed_data_pos) = use_sz;
            compressed_data_pos += sizeof(unsigned char);
            // leave blank for quantizer size
            unsigned char * quantizer_length = compressed_data_pos;
            compressed_data_pos += sizeof(size_t);
            // quantize level elements
            const int quant_radius = 32768;
            // TODO: this is only for 3d
            double C2 = 1 + 3*sqrt(3)/4;
            // double level_eb = eb / (C2 * (target_level + 1));
            double c = sqrt(8);
            fstream file;
            file.open("config");
            if(file){
                file >> c;
                file.close();
            }
            // geometric
            double cc = (1 - c) / (1 - pow(c, target_level + 1));
            double level_eb = cc * eb / C2;
            // arithmetic
            // double cc = 2.0 / (2 + target_level*c) / (target_level + 1);
            // double level_eb = cc * eb / C2;
            // double difference = level_eb * c;
            vector<vector<size_t>> level_dims = init_levels(dims, target_level);
            vector<int> quant_inds(num_elements - num_nodal_elements);
            int quant_count = 0;
            if(use_sz){
                T * data_buffer_pos = data_buffer;
                const T * data_pos = data;
                for(int i=0; i<current_dims[0]; i++){
                    const T * row_data_pos = data_pos;
                    T * row_data_buffer_pos = data_buffer_pos;
                    for(int j=0; j<current_dims[1]; j++){
                        memcpy(row_data_buffer_pos, row_data_pos, current_dims[2] * sizeof(T));
                        row_data_pos += dims[2];
                        row_data_buffer_pos += current_dims[2];
                    }
                    data_pos += dims[1] * dims[2];
                    data_buffer_pos += current_dims[1] * current_dims[2];
                }
                // sz compressing
                size_t sz_compressed_size = 0;
                auto sz_compressed = sz_compress_3d(data_buffer, current_dims[0], current_dims[1], current_dims[2], level_eb, sz_compressed_size, 6);
                // record sz compressed data
                *reinterpret_cast<size_t*>(compressed_data_pos) = sz_compressed_size;
                compressed_data_pos += sizeof(size_t);
                memcpy(compressed_data_pos, sz_compressed, sz_compressed_size);
                compressed_data_pos += sz_compressed_size;
                free(sz_compressed);
            }
            else{
                // resize quantization number
                quant_inds.resize(num_elements);
                vector<size_t> dummy_dims(dims.size(), 0);
                quant_count += quantize_level(data, dims, dummy_dims, level_dims[0], level_eb, quant_radius, quant_inds, quant_count, compressed_data_pos);
            }
            for(int l=1; l<=target_level; l++){
                // geometric
                level_eb *= c;
                // arithmetic
                // level_eb += difference;
                quant_count += quantize_level(data, dims, level_dims[l-1], level_dims[l], level_eb, quant_radius, quant_inds, quant_count, compressed_data_pos);
            }
            // record length for all quantizers
            *reinterpret_cast<size_t*>(quantizer_length) = compressed_data_pos - quantizer_length - sizeof(size_t);
            // encode
    		auto encoder = SZ::HuffmanEncoder<int>();
    		encoder.preprocess_encode(quant_inds, 4*quant_radius);
    		encoder.save(compressed_data_pos);
    		encoder.encode(quant_inds, compressed_data_pos);
    		encoder.postprocess_encode();
        }
		size_t compressed_length = compressed_data_pos - compressed;
        unsigned char * lossless_compressed = NULL;
        size_t lossless_length = sz_lossless_compress(ZSTD_COMPRESSOR, 3, compressed, compressed_length, &lossless_compressed);
        free(compressed);
		compressed_size = lossless_length;
		return lossless_compressed;
	}

	void init(const vector<size_t>& dims){
		size_t buffer_size = default_batch_size * (*max_element(dims.begin(), dims.end())) * sizeof(T);
		// cerr << "buffer_size = " << buffer_size << endl;
		if(data_buffer) free(data_buffer);
		if(correction_buffer) free(correction_buffer);
		if(load_v_buffer) free(load_v_buffer);
		data_buffer = (T *) malloc(data_buffer_size);
		correction_buffer = (T *) malloc(buffer_size);
		load_v_buffer = (T *)malloc(buffer_size);
	}
	// reorder the data to put all the coefficient to the back
	void data_reorder_1D(const T * data_pos, size_t n_nodal, size_t n_coeff, T * nodal_buffer, T * coeff_buffer){
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
	void compute_interpolant_difference_1D(size_t n_coeff, const T * nodal_buffer, T * coeff_buffer){
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
	void decompose_level_1D(T * data_pos, size_t n, T h, bool nodal_row=true){
		size_t n_nodal = (n >> 1) + 1;
		size_t n_coeff = n - n_nodal;
		T * nodal_buffer = data_buffer;
		T * coeff_buffer = data_buffer + n_nodal;
		data_reorder_1D(data_pos, n_nodal, n_coeff, nodal_buffer, coeff_buffer);
		compute_interpolant_difference_1D(n_coeff, nodal_buffer, coeff_buffer);
		if(nodal_row) compute_load_vector_nodal_row(load_v_buffer, n_nodal, n_coeff, h, coeff_buffer);
        else compute_load_vector_coeff_row(load_v_buffer, n_nodal, n_coeff, h, nodal_buffer, coeff_buffer);
		compute_correction(correction_buffer, n_nodal, h, load_v_buffer);
		add_correction(n_nodal, nodal_buffer);
		memcpy(data_pos, data_buffer, n*sizeof(T));
	}
	/*
		2D decomposition
	*/
	// reorder the data to put all the coefficient to the back
	/*
		oxoxo		oooxx		oooxx
		xxxxx	(1)	xxxxx	(2)	oooxx
		oxoxo	=>	oooxx	=>	oooxx
		xxxxx		xxxxx		xxxxx
		oxoxo		oooxx		xxxxx
	*/
	void data_reorder_2D(T * data_pos, size_t n1, size_t n2, size_t stride){
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		size_t n2_nodal = (n2 >> 1) + 1;
		size_t n2_coeff = n2 - n2_nodal;
		T * cur_data_pos = data_pos;
		T * nodal_pos = data_buffer;
		T * coeff_pos = data_buffer + n2_nodal;
		// do reorder (1)
		for(int i=0; i<n1; i++){
			data_reorder_1D(cur_data_pos, n2_nodal, n2_coeff, nodal_pos, coeff_pos);
			memcpy(cur_data_pos, data_buffer, n2 * sizeof(T));
			cur_data_pos += stride;
		}
		if(!(n1 & 1)){
			// n1 is even, change the last coeff row into nodal row
			cur_data_pos -= stride;
			for(int j=0; j<n2; j++){
				cur_data_pos[j] = 2 * cur_data_pos[j] - cur_data_pos[-stride + j];
			}
		}
		// do reorder (2)
		// TODO: change to online processing for memory saving
		switch_rows_2D_by_buffer(data_pos, data_buffer, n1, n2, stride);
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
	void compute_interpolant_difference_2D(T * data_pos, size_t n1, size_t n2, size_t stride){
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		size_t n2_nodal = (n2 >> 1) + 1;
		size_t n2_coeff = n2 - n2_nodal;
		// compute horizontal difference
		const T * nodal_pos = data_pos;
		T * coeff_pos = data_pos + n2_nodal;
		for(int i=0; i<n1_nodal; i++){
			compute_interpolant_difference_1D(n2_coeff, nodal_pos, coeff_pos);
			nodal_pos += stride, coeff_pos += stride;
		}
		// compute vertical difference
		compute_interpolant_difference_2D_vertical(data_pos, n1, n2, stride);
	}	
	// decompose n1 x n2 data into coarse level (n1/2 x n2/2)
	void decompose_level_2D(T * data_pos, size_t n1, size_t n2, T h, size_t stride){
		// cerr << "decompose, h = " << h << endl; 
        size_t n1_nodal = (n1 >> 1) + 1;
        size_t n1_coeff = n1 - n1_nodal;
        size_t n2_nodal = (n2 >> 1) + 1;
        size_t n2_coeff = n2 - n2_nodal;
		data_reorder_2D(data_pos, n1, n2, stride);
		compute_interpolant_difference_2D(data_pos, n1, n2, stride);
        vector<T> w1(n1_nodal);
        vector<T> b1(n1_nodal);
        vector<T> w2(n2_nodal);
        vector<T> b2(n2_nodal);
        precompute_w_and_b(w1.data(), b1.data(), n1_nodal);
        precompute_w_and_b(w2.data(), b2.data(), n2_nodal);
        compute_correction_2D(data_pos, data_buffer, load_v_buffer, n1, n2, n1_nodal, h, stride, w1.data(), b1.data(), w2.data(), b2.data(), default_batch_size);
        apply_correction_batched(data_pos, data_buffer, n1_nodal, stride, n2_nodal, true);
	}
	/*
		2D reorder + vertical reorder
	*/
	void data_reorder_3D(T * data_pos, size_t n1, size_t n2, size_t n3, size_t dim0_stride, size_t dim1_stride){
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		size_t n2_nodal = (n2 >> 1) + 1;
		size_t n2_coeff = n2 - n2_nodal;
		size_t n3_nodal = (n3 >> 1) + 1;
		size_t n3_coeff = n3 - n3_nodal;
		T * cur_data_pos = data_pos;
		// do 2D reorder
		for(int i=0; i<n1; i++){
			data_reorder_2D(cur_data_pos, n2, n3, dim1_stride);
			cur_data_pos += dim0_stride;
		}
		if(!(n1 & 1)){
			// n1 is even, change the last coeff plane into nodal plane
			cur_data_pos -= dim0_stride;
			for(int j=0; j<n2; j++){
				for(int k=0; k<n3; k++){
					cur_data_pos[k] = 2 * cur_data_pos[k] - cur_data_pos[- dim0_stride + k];
				}
				cur_data_pos += dim1_stride;
			}
		}
		cur_data_pos = data_pos;
		// reorder vertically
		for(int j=0; j<n2; j++){
			switch_rows_2D_by_buffer(cur_data_pos, data_buffer, n1, n3, dim0_stride);
			cur_data_pos += dim1_stride;
		}
	}
	/*
		2D computation + vertical computation for coefficient plane 
	*/
	void compute_interpolant_difference_3D(T * data_pos, size_t n1, size_t n2, size_t n3, size_t dim0_stride, size_t dim1_stride){
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		size_t n2_nodal = (n2 >> 1) + 1;
		size_t n2_coeff = n2 - n2_nodal;
		size_t n3_nodal = (n3 >> 1) + 1;
		size_t n3_coeff = n3 - n3_nodal;
		bool even_n2 = (!(n2 & 1));
		bool even_n3 = (!(n3 & 1));
		T * cur_data_pos = data_pos;
 		for(int i=0; i<n1_nodal; i++){
 			compute_interpolant_difference_2D(cur_data_pos, n2, n3, dim1_stride);
 			cur_data_pos += dim0_stride;
 		}
 		// compute vertically
 		const T * nodal_pos = data_pos;
 		T * coeff_pos = data_pos + n1_nodal * dim0_stride;
		for(int i=0; i<n1_coeff; i++){
			// iterate throught coefficient planes along n1
			/*
				data in the coefficient plane
				xxxxx		xxx						xx
				xxxxx		xxx	coeff_nodal_nonal	xx 	coeff_nodal_coeff
				xxxxx	=>	xxx						xx
				xxxxx
				xxxxx		xxx	coeff_coeff_nodal	xx 	coeff_coeff_coeff
							xxx						xx
			*/
            const T * nodal_nodal_nodal_pos = nodal_pos;
            T * coeff_nodal_nodal_pos = coeff_pos;
            T * coeff_nodal_coeff_pos = coeff_pos + n3_nodal;
            T * coeff_coeff_nodal_pos = coeff_pos + n2_nodal * dim1_stride;
            T * coeff_coeff_coeff_pos = coeff_coeff_nodal_pos + n3_nodal;
            // TODO: optimize average computation
            for(int j=0; j<n2_coeff; j++){
            	for(int k=0; k<n3_coeff; k++){
	                // coeff_nodal_nonal
	                coeff_nodal_nodal_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k]) / 2;
	                // coeff_nodal_coeff
	                coeff_nodal_coeff_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k] +
	                								nodal_nodal_nodal_pos[k + 1] + nodal_nodal_nodal_pos[dim0_stride + k + 1]) / 4;
	                // coeff_coeff_nodal
	                coeff_coeff_nodal_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k] +
	                								nodal_nodal_nodal_pos[k + dim1_stride] + nodal_nodal_nodal_pos[dim0_stride + k + dim1_stride]) / 4;
	                // coeff_coeff_coeff
	                coeff_coeff_coeff_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k] +
	                								nodal_nodal_nodal_pos[k + 1] + nodal_nodal_nodal_pos[dim0_stride + k + 1] +
	                								nodal_nodal_nodal_pos[k + dim1_stride] + nodal_nodal_nodal_pos[dim0_stride + k + dim1_stride] + 
	                								nodal_nodal_nodal_pos[k + dim1_stride + 1] + nodal_nodal_nodal_pos[dim0_stride + k + dim1_stride + 1]) / 8;
                }
	            // compute the last (or second last if n3 is even) coeff_*_nodal column
	            coeff_nodal_nodal_pos[n3_coeff] -= (nodal_nodal_nodal_pos[n3_coeff] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff]) / 2;
                coeff_coeff_nodal_pos[n3_coeff] -= (nodal_nodal_nodal_pos[n3_coeff] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff] +
                								nodal_nodal_nodal_pos[n3_coeff + dim1_stride] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + dim1_stride]) / 4;
	            if(even_n3){
	            	// compute the last coeff_*_nodal column if n3 is even
		            coeff_nodal_nodal_pos[n3_coeff + 1] -= (nodal_nodal_nodal_pos[n3_coeff + 1] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]) / 2;
	                coeff_coeff_nodal_pos[n3_coeff + 1] -= (nodal_nodal_nodal_pos[n3_coeff + 1] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1] +
	                								nodal_nodal_nodal_pos[n3_coeff + 1 + dim1_stride] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1 + dim1_stride]) / 4;
	            }
	            coeff_nodal_nodal_pos += dim1_stride;
	            coeff_nodal_coeff_pos += dim1_stride;
	            coeff_coeff_nodal_pos += dim1_stride;
	            coeff_coeff_coeff_pos += dim1_stride;
	            nodal_nodal_nodal_pos += dim1_stride;
            }
            // compute the last (or second last if n2 is even) coeff_nodal_* row
            {
            	for(int k=0; k<n3_coeff; k++){
            		coeff_nodal_nodal_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k]) / 2;
	                coeff_nodal_coeff_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k] +
	                								nodal_nodal_nodal_pos[k + 1] + nodal_nodal_nodal_pos[dim0_stride + k + 1]) / 4;
            	}
        		coeff_nodal_nodal_pos[n3_coeff] -= (nodal_nodal_nodal_pos[n3_coeff] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff]) / 2;
            	if(even_n3){
            		coeff_nodal_nodal_pos[n3_coeff + 1] -= (nodal_nodal_nodal_pos[n3_coeff + 1] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]) / 2;
            	}
                coeff_nodal_nodal_pos += dim1_stride;
                coeff_nodal_coeff_pos += dim1_stride;
                coeff_coeff_nodal_pos += dim1_stride;
                coeff_coeff_coeff_pos += dim1_stride;
                nodal_nodal_nodal_pos += dim1_stride;
        	}
            if(even_n2){
                // compute the last coeff_nodal_* row if n2 is even
            	for(int k=0; k<n3_coeff; k++){
            		coeff_nodal_nodal_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k]) / 2;
	                coeff_nodal_coeff_pos[k] -= (nodal_nodal_nodal_pos[k] + nodal_nodal_nodal_pos[dim0_stride + k] +
	                								nodal_nodal_nodal_pos[k + 1] + nodal_nodal_nodal_pos[dim0_stride + k + 1]) / 4;
            	}
        		coeff_nodal_nodal_pos[n3_coeff] -= (nodal_nodal_nodal_pos[n3_coeff] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff]) / 2;
            	if(even_n3){
            		coeff_nodal_nodal_pos[n3_coeff + 1] -= (nodal_nodal_nodal_pos[n3_coeff + 1] + nodal_nodal_nodal_pos[dim0_stride + n3_coeff + 1]) / 2;
            	}
            }
            nodal_pos += dim0_stride;
            coeff_pos += dim0_stride;
		}
	}

	// decompse n1 x n2 x n3 data into coarse level (n1/2 x n2/2 x n3/2)
	void decompose_level_3D(T * data_pos, size_t n1, size_t n2, size_t n3, T h, size_t dim0_stride, size_t dim1_stride){
		data_reorder_3D(data_pos, n1, n2, n3, dim0_stride, dim1_stride);
		compute_interpolant_difference_3D(data_pos, n1, n2, n3, dim0_stride, dim1_stride);
        size_t n1_nodal = (n1 >> 1) + 1;
        size_t n2_nodal = (n2 >> 1) + 1;
        size_t n3_nodal = (n3 >> 1) + 1;
        compute_correction_3D(data_pos, data_buffer, load_v_buffer, n1, n2, n3, n1_nodal, h, dim0_stride, dim1_stride, default_batch_size);
        T * nodal_pos = data_pos;
        const T * correction_pos = data_buffer;
        for(int i=0; i<n1_nodal; i++){
            apply_correction_batched(nodal_pos, correction_pos, n2_nodal, dim1_stride, n3_nodal, true);
            nodal_pos += dim0_stride;
            correction_pos += n2_nodal * n3_nodal;
        }
	}
};


}

#endif