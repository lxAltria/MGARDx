#ifndef _MGARD_ADAPTIVE_RESOLUTION_HPP
#define _MGARD_ADAPTIVE_RESOLUTION_HPP

#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include "decompose.hpp"
#include "recompose.hpp"

namespace MGARD{

using namespace std;

template <class T>
class AR_Decomposer{
public:
	AR_Decomposer(const vector<size_t>& dims, uint block_size, bool method=0) 
		: method(method), block_size(block_size), dims(dims) {
		strides = vector<size_t>(dims.size());
		size_t stride = 1;
		for(int i=dims.size()-1; i>=0; i--){
			strides[i] = stride;
			stride *= dims[i];
		} 
		num_blocks = vector<size_t>(dims.size());
		for(int i=0; i<dims.size(); i++){
			num_blocks[i] = (dims[i] - 1) / block_size + 1;
		}
		// simplify method index to option
    	option = (method == 0);
	}

	// 3D decomposition
	void decompose(T * data){
		auto nx = num_blocks[0];
		auto ny = num_blocks[1];
		auto nz = num_blocks[2];
	    T * x_data_pos = data;
	    vector<size_t> block_dims(3, 0);
	    for(int i=0; i<nx; i++){
	    	T * y_data_pos = x_data_pos;
	    	block_dims[0] = (i < nx - 1) ? block_size : (dims[0] - i*block_size);
	    	for(int j=0; j<ny; j++){
	    		T * z_data_pos = y_data_pos;
		    	block_dims[1] = (j < ny - 1) ? block_size : (dims[1] - j*block_size);
	    		for(int k=0; k<nz; k++){
			    	block_dims[2] = (k < nz - 1) ? block_size : (dims[2] - k*block_size);
	    			cout << "dimensions = " << block_dims[0] << " " << block_dims[1] << " " << block_dims[2] << endl;
				    MGARD::Decomposer<T> decomposer;
			    	decomposer.decompose(z_data_pos, block_dims, target_level, option, strides);
			    	z_data_pos += block_size;
	    		}
	    		y_data_pos += block_size * dims[2];
	    	}
	    	x_data_pos += block_size * dims[1] * dims[2];
	    }		
	}

	void progressive_recompose(T * data){
		auto nx = num_blocks[0];
		auto ny = num_blocks[1];
		auto nz = num_blocks[2];
	    T * x_data_pos = data;
	    vector<size_t> block_dims(3, 0);
	    int block_id = 0;
	    for(int i=0; i<nx; i++){
	    	T * y_data_pos = x_data_pos;
	    	block_dims[0] = (i < nx - 1) ? block_size : (dims[0] - i*block_size);
	    	for(int j=0; j<ny; j++){
	    		T * z_data_pos = y_data_pos;
		    	block_dims[1] = (j < ny - 1) ? block_size : (dims[1] - j*block_size);
	    		for(int k=0; k<nz; k++){
			    	block_dims[2] = (k < nz - 1) ? block_size : (dims[2] - k*block_size);
			    	if(converged[block_id]) continue;
				    MGARD::Recomposer<T> recomposer;
				    init_levels(block_dims, max_level - current_level[block_id]);
				    auto next_dims = level_dims[1];
				    int target_level = 1;
				    // check maximal level for boundary?
			    	recomposer.recompose(z_data_pos, block_dims, target_level, option, strides);
			    	current_level[block_id] ++;
			    	// placeholder for testing
			    	analysis_result[block_id].push_back(0);
			    	// analysis_result[block_id].push_back(perform_analysis(z_data_pos, next_dims, strides));
			    	// converged[block_id] = check_convergence(analysis_result[block_id]);
			    	block_id ++;
			    	z_data_pos += block_size;
	    		}
	    		y_data_pos += block_size * dims[2];
	    	}
	    	x_data_pos += block_size * dims[1] * dims[2];
	    }
	}

private:
	uint method = 0; // 0: hierarchical; 1: orthogonal (MGARD)
	bool option = false;
	uint block_size = 0;
	vector<size_t> dims;
	vector<size_t> num_blocks;
	vector<size_t> strides;
	Decomposer<T> decomposer;
	Recomposer<T> recomposer;
	vector<uint> current_level;
	vector<bool> converged;
	vector<vector<double>> analysis_result;
}
#endif