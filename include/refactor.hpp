#ifndef _REFACTOR_HPP
#define _REFACTOR_HPP

#include <vector>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "utils.hpp"

namespace MGARD{

using namespace std;

// interleave coeffcients in the finer level
/*
@params data: decomposed data
@params dims: original dims (to compute offset)
@params dims_fine: fine level dimensions
@params dims_coarse: coarse level dimensions
@params buffer: data buffer for level coefficients
@params max_err: max error in current level
*/
template <class T>
void interleave_level_coefficients_3d(const T * data, const vector<size_t>& dims, const vector<size_t>& dims_fine, const vector<size_t>& dims_coasre, T * buffer, T& max_err){
    size_t dim0_offset = dims[1] * dims[2];
    size_t dim1_offset = dims[2];
    max_err = 0;
    size_t count = 0;
    for(int i=0; i<dims_fine[0]; i++){
        for(int j=0; j<dims_fine[1]; j++){
            for(int k=0; k<dims_fine[2]; k++){
                if((i < dims_coasre[0]) && (j < dims_coasre[1]) && (k < dims_coasre[2]))
                    continue;
                buffer[count] = data[i*dim0_offset + j*dim1_offset + k];
                T err = fabs(buffer[count]);
                if(err > max_err) max_err = err;
                count ++;
            }
        }
    }
}
template <class T>
void interleave_level_coefficients(const T * data, const vector<size_t>& dims, const vector<size_t>& dims_fine, const vector<size_t>& dims_coasre, T * buffer, T& max_err){
    switch(dims.size()){
        case 3:
            interleave_level_coefficients_3d(data, dims, dims_fine, dims_coasre, buffer, max_err);
            break;
        default:
            cerr << "Other dimensions are not supported" << endl;
    }
}
// refactor level-centric decomposed data in hierarchical fashion
/*
@params data: decomposed data
@params target_level: decomposed level
@params dims: data dimensions
*/
template <class T>
vector<unsigned char*> level_centric_data_refactor(const T * data, int target_level, const vector<size_t>& dims){
    int max_level = log2(*min_element(dims.begin(), dims.end()));
    if(target_level > max_level) target_level = max_level;
    auto level_dims = init_levels(dims, target_level);
    unsigned char * metadata = (unsigned char *) malloc((target_level + 1) * (sizeof(size_t) + sizeof(T)));
    size_t * level_elements = reinterpret_cast<size_t*>(metadata);
    T * level_error_bounds = reinterpret_cast<T*>(metadata + (target_level + 1) * sizeof(size_t));
    level_elements[0] = 1;
    for(int j=0; j<dims.size(); j++){
        level_elements[0] *= level_dims[0][j];
    }
    size_t pre_num_elements = level_elements[0];
    for(int i=1; i<=target_level; i++){
        size_t num_elements = 1;
        for(int j=0; j<dims.size(); j++){
            num_elements *= level_dims[i][j];
        }
        level_elements[i] = num_elements - pre_num_elements;
        pre_num_elements = num_elements;
    }
    vector<unsigned char *> level_components;
    level_components.push_back(metadata);
    vector<size_t> dims_dummy(dims.size(), 0);
    unsigned char * buffer = (unsigned char *) malloc(level_elements[0] * sizeof(T));
    interleave_level_coefficients(data, dims, level_dims[0], dims_dummy, reinterpret_cast<T*>(buffer), level_error_bounds[0]);
    level_components.push_back(buffer);
    fflush(stdout);
    for(int i=1; i<=target_level; i++){
        unsigned char * buffer = (unsigned char *) malloc(level_elements[i] * sizeof(T));
        interleave_level_coefficients(data, dims, level_dims[i], level_dims[i - 1], reinterpret_cast<T*>(buffer), level_error_bounds[i]);
        level_components.push_back(buffer);
    }
    return level_components;
}

// put coeffcients in the finer level to the correct position
/*
@params buffer: data buffer for level coefficients
@params dims: original dims (to compute offset)
@params dims_fine: fine level dimensions
@params dims_coarse: coarse level dimensions
@params data: decomposed data
*/
template <class T>
void reposition_level_coefficients_3d(const T * buffer, const vector<size_t>& dims, const vector<size_t>& dims_fine, const vector<size_t>& dims_coasre, T * data){
    size_t dim0_offset = dims[1] * dims[2];
    size_t dim1_offset = dims[2];
    int count = 0;
    for(int i=0; i<dims_fine[0]; i++){
        for(int j=0; j<dims_fine[1]; j++){
            for(int k=0; k<dims_fine[2]; k++){
                if((i < dims_coasre[0]) && (j < dims_coasre[1]) && (k < dims_coasre[2]))
                    continue;
                data[i*dim0_offset + j*dim1_offset + k] = buffer[count ++];
            }
        }
    }
}
template <class T>
void reposition_level_coefficients(const T * buffer, const vector<size_t>& dims, const vector<size_t>& dims_fine, const vector<size_t>& dims_coasre, T * data){
    switch(dims.size()){
        case 3:
            reposition_level_coefficients_3d(buffer, dims, dims_fine, dims_coasre, data);
            break;
        default:
            cerr << "Other dimensions are not supported" << endl;
    }
}
// reposition level-centric decomposed data for recomposition
/*
@params data: decomposed data
@params target_level: decomposed level
@params target_recompose_level: recomposed level
@params dims: data dimensions, modified after calling the function
*/
template <class T>
T * level_centric_data_reposition(const vector<unsigned char*>& level_components, int target_level, int target_recompose_level, vector<size_t>& dims){
    auto level_dims = init_levels(dims, target_level);
    unsigned char * metadata = level_components[0];
    size_t * level_elements = reinterpret_cast<size_t*>(metadata);
    T * level_error_bounds = reinterpret_cast<T*>(metadata + (target_level + 1) * sizeof(size_t));    
    size_t num_elements = 1;
    for(int j=0; j<dims.size(); j++){
        dims[j] = level_dims[target_recompose_level][j];
        num_elements *= dims[j];
    }
    T * data = (T *) malloc(num_elements * sizeof(T));
    vector<size_t> dims_dummy(dims.size(), 0);
    reposition_level_coefficients(reinterpret_cast<T*>(level_components[1]), dims, level_dims[0], dims_dummy, data);
    for(int i=1; i<=target_recompose_level; i++){
        reposition_level_coefficients(reinterpret_cast<T*>(level_components[i + 1]), dims, level_dims[i], level_dims[i - 1], data);        
    }
    return data;
}

}

#endif