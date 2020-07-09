#ifndef _REFACTOR_HPP
#define _REFACTOR_HPP

#include <vector>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <bitset>
#include <iomanip>
#include "utils.hpp"
#include "data_org.hpp"
#include "data_enc.hpp"
#include "error_est.hpp"

namespace REFACTOR{

using namespace std;
using namespace MGARD;

// size of segment: default 4 MB
const int seg_size = 4;

// metadata for refactored data
template <class T>
class Metadata{
public:
    Metadata(){}
    Metadata(int target_level){
        level_elements = vector<size_t>(target_level + 1);
        level_error_bounds = vector<T>(target_level + 1);
    }
    vector<size_t> level_elements;    // number of elements in each multigrid level
    vector<T> level_error_bounds;     // max errors in each multigrid level
    vector<vector<size_t>> component_sizes; // size of each component in each level
    vector<vector<double>> max_e;         // max_e[i][j]: max error of level i using the first (j+1) bit-planes
    vector<vector<double>> mse;         // mse[i][j]: mse of level i using the first (j+1) bit-planes
    vector<vector<unsigned char>> bitplane_indictors;   // indicator for bitplane encoding
    vector<int> order;                  // order of bitplane placement
    bool max_e_estimator = false;
    bool mse_estimator = false;
    int option = ENCODING_DEFAULT;
    int encoded_bitplanes = 32;

    void init_encoded_sizes(){
        component_sizes.clear();
        for(int i=0; i<level_elements.size(); i++){
            component_sizes.push_back(vector<size_t>());
        }
    }
    void init_bitplane_indicators(){
        bitplane_indictors.clear();
        for(int i=0; i<level_elements.size(); i++){
            bitplane_indictors.push_back(vector<unsigned char>());
        }        
    }
    size_t size(){
        size_t metadata_size = 
                sizeof(int) + sizeof(int) + sizeof(size_t)  // option + encoded_bitplanes + number of levels
                + level_elements.size() * sizeof(size_t)    // level_elements
                + level_error_bounds.size() * sizeof(T)     // level_eb
                + sizeof(unsigned char) + sizeof(unsigned char) // estimator flags
                + sizeof(size_t) + order.size() * sizeof(int) // order
                ;
        for(int i=0; i<component_sizes.size(); i++){
            metadata_size += sizeof(size_t) + component_sizes[i].size() * sizeof(size_t);
        }
        for(int i=0; i<bitplane_indictors.size(); i++){
            metadata_size += sizeof(size_t) + bitplane_indictors[i].size() * sizeof(unsigned char);
        }
        if(max_e_estimator) {
            for(int i=0; i<max_e.size(); i++){
                metadata_size += sizeof(size_t) + max_e[i].size() * sizeof(double);
            }
        }
        if(mse_estimator) {
            for(int i=0; i<mse.size(); i++){
                metadata_size += sizeof(size_t) + mse[i].size() * sizeof(double);
            }
        }
        return metadata_size;
    }
    template <class T1>
    size_t serialize_level_vectors(const vector<vector<T1>>& level_vecs, unsigned char * buffer_pos){
        unsigned char const * const start = buffer_pos;
        for(int i=0; i<level_vecs.size(); i++){
            *reinterpret_cast<size_t*>(buffer_pos) = level_vecs[i].size();
            buffer_pos += sizeof(size_t);
            memcpy(buffer_pos, level_vecs[i].data(), level_vecs[i].size() * sizeof(T1));
            buffer_pos += level_vecs[i].size() * sizeof(T1);
        }
        return buffer_pos - start;
    }
    unsigned char * serialize(){
        unsigned char * buffer = (unsigned char *) malloc(size());
        unsigned char * buffer_pos = buffer;
        *reinterpret_cast<int*>(buffer_pos) = option;
        buffer_pos += sizeof(int);
        *reinterpret_cast<int*>(buffer_pos) = encoded_bitplanes;
        buffer_pos += sizeof(int);
        size_t num_levels = level_elements.size();
        *reinterpret_cast<size_t*>(buffer_pos) = num_levels;
        buffer_pos += sizeof(size_t);
        memcpy(buffer_pos, level_elements.data(), num_levels * sizeof(size_t));
        buffer_pos += num_levels * sizeof(size_t);
        memcpy(buffer_pos, level_error_bounds.data(), num_levels * sizeof(T));
        buffer_pos += num_levels * sizeof(T);
        *reinterpret_cast<size_t*>(buffer_pos) = order.size();
        buffer_pos += sizeof(size_t);
        memcpy(buffer_pos, order.data(), order.size() * sizeof(int));
        buffer_pos += order.size() * sizeof(int);
        *buffer_pos = mse_estimator;
        buffer_pos += sizeof(unsigned char);
        *buffer_pos = max_e_estimator;
        buffer_pos += sizeof(unsigned char);
        buffer_pos += serialize_level_vectors(component_sizes, buffer_pos);
        buffer_pos += serialize_level_vectors(bitplane_indictors, buffer_pos);
        if(max_e_estimator){
            buffer_pos += serialize_level_vectors(max_e, buffer_pos);
        }
        if(mse_estimator){
            buffer_pos += serialize_level_vectors(mse, buffer_pos);
        }
        return buffer;
    }
    // auto increment buffer_pos
    template <class T1>
    vector<vector<T1>> deserialize_level_vectors(const unsigned char *& buffer_pos, size_t num_levels){
        vector<vector<T1>> level_vecs;
        for(int i=0; i<num_levels; i++){
            size_t num = *reinterpret_cast<const size_t*>(buffer_pos);
            buffer_pos += sizeof(size_t);
            vector<T1> level_vec = vector<T1>(reinterpret_cast<const T1 *>(buffer_pos), reinterpret_cast<const T1 *>(buffer_pos) + num);
            level_vecs.push_back(level_vec);
            buffer_pos += num * sizeof(T1);
        }
        return level_vecs;
    }
    void deserialize(const unsigned char * serialized_data){
        const unsigned char * buffer_pos = serialized_data;
        option = *reinterpret_cast<const int*>(buffer_pos);
        buffer_pos += sizeof(int);
        encoded_bitplanes = *reinterpret_cast<const int*>(buffer_pos);
        buffer_pos += sizeof(int);
        size_t num_levels = *reinterpret_cast<const size_t*>(buffer_pos);
        buffer_pos += sizeof(size_t);
        level_elements = vector<size_t>(reinterpret_cast<const size_t *>(buffer_pos), reinterpret_cast<const size_t *>(buffer_pos) + num_levels);
        buffer_pos += num_levels * sizeof(size_t);
        level_error_bounds = vector<T>(reinterpret_cast<const T *>(buffer_pos), reinterpret_cast<const T *>(buffer_pos) + num_levels);
        buffer_pos += num_levels * sizeof(T);
        size_t order_size = *reinterpret_cast<const size_t*>(buffer_pos);
        buffer_pos += sizeof(size_t);
        order = vector<int>(reinterpret_cast<const int *>(buffer_pos), reinterpret_cast<const int *>(buffer_pos) + order_size);
        buffer_pos += order_size * sizeof(int);
        mse_estimator = *buffer_pos;
        buffer_pos += sizeof(unsigned char);
        max_e_estimator = *buffer_pos;
        buffer_pos += sizeof(unsigned char);
        // deserialize_level_vectors has auto increment for buffer_pos
        component_sizes = deserialize_level_vectors<size_t>(buffer_pos, num_levels);
        bitplane_indictors = deserialize_level_vectors<unsigned char>(buffer_pos, num_levels);
        if(max_e_estimator){
            // deserialize_level_vectors has auto increment for buffer_pos
            max_e = deserialize_level_vectors<double>(buffer_pos, num_levels);
        }
        if(mse_estimator){
            // deserialize_level_vectors has auto increment for buffer_pos
            mse = deserialize_level_vectors<double>(buffer_pos, num_levels);
        }
    }
    void to_file(const string& filename){
        auto serialized_data = serialize();
        writefile(filename.c_str(), serialized_data, size());
        free(serialized_data);
    }
    void from_file(const string& filename){
        size_t num_bytes = 0;
        unsigned char * buffer = readfile_pointer<unsigned char>(filename.c_str(), num_bytes);
        deserialize(buffer);
        free(buffer);
    }
};

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
void interleave_level_coefficients(const T * data, const vector<size_t>& dims, const vector<size_t>& dims_fine, const vector<size_t>& dims_coasre, T * buffer){
    switch(dims.size()){
        case 3:
            interleave_level_coefficients_3d(data, dims, dims_fine, dims_coasre, buffer);
            break;
        default:
            cerr << "Other dimensions are not supported" << endl;
    }
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
void reposition_level_coefficients(const T * buffer, const vector<size_t>& dims, const vector<size_t>& dims_fine, const vector<size_t>& dims_coasre, T * data){
    switch(dims.size()){
        case 3:
            reposition_level_coefficients_3d(buffer, dims, dims_fine, dims_coasre, data);
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
@params metadata: metadata object
*/
template <class T>
vector<vector<unsigned char*>> level_centric_data_refactor(const T * data, int target_level, const vector<size_t>& dims, Metadata<T>& metadata){
    int max_level = log2(*min_element(dims.begin(), dims.end()));
    if(target_level > max_level) target_level = max_level;
    auto level_dims = init_levels(dims, target_level);
    vector<size_t>& level_elements = metadata.level_elements;
    vector<T>& level_error_bounds = metadata.level_error_bounds;
    // compute level elements
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
    // init bitplane indicators
    metadata.init_bitplane_indicators();
    // init metadata sizes recorder
    metadata.init_encoded_sizes();
    metadata.max_e_estimator = true;
    // turn on mse
    metadata.mse_estimator = true;
    // record all level components
    vector<vector<unsigned char*>> level_components;
    vector<size_t> dims_dummy(dims.size(), 0);
    for(int i=0; i<=target_level; i++){
        // cout << "encoding level " << target_level - i << endl;
        const vector<size_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
        unsigned char * buffer = (unsigned char *) malloc(level_elements[i] * sizeof(T));
        // extract components for each level
        interleave_level_coefficients(data, dims, level_dims[i], prev_dims, reinterpret_cast<T*>(buffer));

        string outfile("decomposed_level_");
        writefile((outfile + to_string(target_level - i) + ".dat").c_str(), reinterpret_cast<T*>(buffer), level_elements[i]);

        level_error_bounds[i] = record_level_max_value(reinterpret_cast<T*>(buffer), level_elements[i]);
        if(level_elements[i] * sizeof(T) < seg_size){
            if(metadata.max_e_estimator){
                auto level_max_e = vector<double>(1, 0);
                metadata.max_e.push_back(level_max_e);
            }
            if(metadata.mse_estimator){
                auto level_mse = vector<double>(1, 0);
                metadata.mse.push_back(level_mse);
            }
            vector<unsigned char*> tiny_level;
            tiny_level.push_back(buffer);
            level_components.push_back(tiny_level);
            metadata.component_sizes[i].push_back(level_elements[i] * sizeof(T));
        }
        else{
            // identify exponent of max element
            int level_exp = 0;
            frexp(level_error_bounds[i], &level_exp);
            cout << "level " << i << " max err = " << level_error_bounds[i] << ", exp = " << level_exp << endl;
            if(metadata.max_e_estimator){
                auto level_max_e = record_level_max_e(reinterpret_cast<T*>(buffer), level_elements[i], metadata.encoded_bitplanes, level_exp);
                metadata.max_e.push_back(level_max_e);
            }
            if(metadata.mse_estimator){
                auto level_mse = record_level_mse(reinterpret_cast<T*>(buffer), level_elements[i], metadata.encoded_bitplanes, level_exp);
                metadata.mse.push_back(level_mse);
            }
            // intra-level progressive encoding
            if(metadata.option == ENCODING_DEFAULT){
                struct timespec start, end;
                int err = clock_gettime(CLOCK_REALTIME, &start);
                auto intra_level_components = progressive_encoding(reinterpret_cast<T*>(buffer), level_elements[i], level_exp, metadata.encoded_bitplanes, metadata.component_sizes[i]);
                level_components.push_back(intra_level_components);
                err = clock_gettime(CLOCK_REALTIME, &end);
                cout << "Byteplane encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
            }
            else if(metadata.option == ENCODING_RLE){
                auto intra_level_components = progressive_encoding_with_rle_compression(reinterpret_cast<T*>(buffer), level_elements[i], level_exp, metadata.encoded_bitplanes, metadata.component_sizes[i]);
                level_components.push_back(intra_level_components);
            }
            else if(metadata.option == ENCODING_HYBRID){
                vector<unsigned char>& bitplane_indictor = metadata.bitplane_indictors[i];
                auto intra_level_components = progressive_hybrid_encoding(reinterpret_cast<T*>(buffer), level_elements[i], level_exp, metadata.encoded_bitplanes, metadata.component_sizes[i], bitplane_indictor);
                level_components.push_back(intra_level_components);
            }
            // release extracted component
            free(buffer);
        }
    }
    return level_components;
}

// reposition level-centric decomposed data for recomposition
/*
@params data: decomposed data
@params metadata: metadata object
@params target_level: decomposed level
@params target_recompose_level: recomposed level
@params dims: data dimensions, modified after calling the function
*/
template <class T>
T * level_centric_data_reposition(const vector<vector<const unsigned char*>>& level_components, const Metadata<T>& metadata, int target_level, int target_recompose_level, const vector<int>& intra_recompose_level, vector<size_t>& dims){
    auto level_dims = init_levels(dims, target_level);
    const vector<size_t>& level_elements = metadata.level_elements;
    const vector<T>& level_error_bounds = metadata.level_error_bounds;
    size_t num_elements = 1;
    for(int j=0; j<dims.size(); j++){
        dims[j] = level_dims[target_recompose_level][j];
        num_elements *= dims[j];
    }
    T * data = (T *) malloc(num_elements * sizeof(T));
    vector<size_t> dims_dummy(dims.size(), 0);
    double total_mse = 0;
    const int mse_factor = 1 << dims.size();
    // reposition_level_coefficients(reinterpret_cast<T*>(level_components[0][0]), dims, level_dims[0], dims_dummy, data);
    for(int i=0; i<=target_recompose_level; i++){
        const vector<size_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
        if(level_elements[i] * sizeof(T) < seg_size){
            reposition_level_coefficients(reinterpret_cast<const T *>(level_components[i][0]), dims, level_dims[i], prev_dims, data);
        }
        else{
            cout << "decoding level " << level_elements.size() - 1 - i << ", size of components = " << level_components[i].size() << endl;
            // identify exponent of max element
            int level_exp = 0;
            frexp(level_error_bounds[i], &level_exp);
            T * buffer = NULL;
            int encoded_bitplanes = intra_recompose_level[i] ? intra_recompose_level[i] : metadata.encoded_bitplanes;
            if(encoded_bitplanes > metadata.encoded_bitplanes) encoded_bitplanes = metadata.encoded_bitplanes;            
            // cout << "encoded_bitplanes = " << encoded_bitplanes << endl;
            if(metadata.mse_estimator){
                cout << "MSE: " << i << " " << encoded_bitplanes - 1 << ":" << metadata.mse[i][encoded_bitplanes - 1] << endl; 
                cout << total_mse << ": " << total_mse * mse_factor << " + " << metadata.mse[i][encoded_bitplanes - 1] << " = ";
                total_mse = total_mse * mse_factor + metadata.mse[i][encoded_bitplanes - 1];
                cout << total_mse << endl;
            }
            // intra-level progressive decoding
            if(metadata.option == ENCODING_DEFAULT){
                // struct timespec start, end;
                // int err = clock_gettime(CLOCK_REALTIME, &start);
                buffer = progressive_decoding<T>(level_components[i], level_elements[i], level_exp, encoded_bitplanes);
                // err = clock_gettime(CLOCK_REALTIME, &end);
                // cout << "Byteplane decoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
            }
            else if(metadata.option == ENCODING_RLE){
                buffer = progressive_decoding_with_rle_compression<T>(level_components[i], level_elements[i], level_exp, encoded_bitplanes);
            }
            else if(metadata.option == ENCODING_HYBRID){
                const vector<unsigned char>& bitplane_indictor = metadata.bitplane_indictors[i];
                buffer = progressive_hybrid_decoding<T>(level_components[i], level_elements[i], level_exp, encoded_bitplanes, bitplane_indictor);
            }
            reposition_level_coefficients(buffer, dims, level_dims[i], prev_dims, data);

            string outfile("reconstructed_level_");
            writefile((outfile + to_string(level_elements.size() - 1 - i) + ".dat").c_str(), reinterpret_cast<const T *>(buffer), level_elements[i]);
            // release reconstructed component
            free(buffer);
        }

    }
    total_mse /= num_elements;
    cout << "total MSE = " << total_mse << endl;
    return data;
}

}

#endif