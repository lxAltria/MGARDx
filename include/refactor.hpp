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

// BitEncoder and BitDecoder are modified from ZFP-0.3.1
// TODO: switch to ZFP-0.5.4 for better performance
/*
** Copyright (c) 2014, Lawrence Livermore National Security, LLC.
** Produced at the Lawrence Livermore National Laboratory.
** Written by Peter Lindstrom.
** LLNL-CODE-663824.
** All rights reserved.
*/
// A class to write data bit by bit
// Currently does not consider the case when encoding size is larger than capacity
#ifdef __GNUC__
  #define align_(n) __attribute__((aligned(n)))
#else
  #define align_(n)
#endif
class BitEncoder{
protected:
    // bit reversal table used by writer
    static unsigned char reverse(unsigned char x){
        static const unsigned char lut[] align_(256) = {
        0x00,0x80,0x40,0xc0,0x20,0xa0,0x60,0xe0,0x10,0x90,0x50,0xd0,0x30,0xb0,0x70,0xf0,
        0x08,0x88,0x48,0xc8,0x28,0xa8,0x68,0xe8,0x18,0x98,0x58,0xd8,0x38,0xb8,0x78,0xf8,
        0x04,0x84,0x44,0xc4,0x24,0xa4,0x64,0xe4,0x14,0x94,0x54,0xd4,0x34,0xb4,0x74,0xf4,
        0x0c,0x8c,0x4c,0xcc,0x2c,0xac,0x6c,0xec,0x1c,0x9c,0x5c,0xdc,0x3c,0xbc,0x7c,0xfc,
        0x02,0x82,0x42,0xc2,0x22,0xa2,0x62,0xe2,0x12,0x92,0x52,0xd2,0x32,0xb2,0x72,0xf2,
        0x0a,0x8a,0x4a,0xca,0x2a,0xaa,0x6a,0xea,0x1a,0x9a,0x5a,0xda,0x3a,0xba,0x7a,0xfa,
        0x06,0x86,0x46,0xc6,0x26,0xa6,0x66,0xe6,0x16,0x96,0x56,0xd6,0x36,0xb6,0x76,0xf6,
        0x0e,0x8e,0x4e,0xce,0x2e,0xae,0x6e,0xee,0x1e,0x9e,0x5e,0xde,0x3e,0xbe,0x7e,0xfe,
        0x01,0x81,0x41,0xc1,0x21,0xa1,0x61,0xe1,0x11,0x91,0x51,0xd1,0x31,0xb1,0x71,0xf1,
        0x09,0x89,0x49,0xc9,0x29,0xa9,0x69,0xe9,0x19,0x99,0x59,0xd9,0x39,0xb9,0x79,0xf9,
        0x05,0x85,0x45,0xc5,0x25,0xa5,0x65,0xe5,0x15,0x95,0x55,0xd5,0x35,0xb5,0x75,0xf5,
        0x0d,0x8d,0x4d,0xcd,0x2d,0xad,0x6d,0xed,0x1d,0x9d,0x5d,0xdd,0x3d,0xbd,0x7d,0xfd,
        0x03,0x83,0x43,0xc3,0x23,0xa3,0x63,0xe3,0x13,0x93,0x53,0xd3,0x33,0xb3,0x73,0xf3,
        0x0b,0x8b,0x4b,0xcb,0x2b,0xab,0x6b,0xeb,0x1b,0x9b,0x5b,0xdb,0x3b,0xbb,0x7b,0xfb,
        0x07,0x87,0x47,0xc7,0x27,0xa7,0x67,0xe7,0x17,0x97,0x57,0xd7,0x37,0xb7,0x77,0xf7,
        0x0f,0x8f,0x4f,0xcf,0x2f,0xaf,0x6f,0xef,0x1f,0x9f,0x5f,0xdf,0x3f,0xbf,0x7f,0xff,
        };
        return lut[x];
    }    
public:
    BitEncoder(unsigned char * array, size_t size){
        start = array;
        current = array;
        capacity = size;
        buffer = 1u;
    }
    void encode(bool b){
        if(current - start > capacity){
            cout << current - start << " " << capacity << "????\n";
            exit(0);
        }
        buffer = (buffer << 1u) + b;
        if(buffer >= 0x100u){
            *current++ = reverse(buffer - 0x100u);
            buffer = 1u;
        }
    }
    void flush(){
        while (buffer != 1u) encode(false);
    }
    size_t size(){ return (buffer == 1u) ? current - start : current - start + 1; }
private:
    unsigned char * start = NULL;
    unsigned char * current = NULL;
    size_t capacity = 0;
    uint buffer = 1u;
};
// A class to read data bit by bit
// Currently does not consider the case when decoding size is larger than capacity
class BitDecoder{
public:
    BitDecoder(unsigned char const * array, size_t size){
        start = array;
        current = array;
        capacity = size;
        buffer = 1u;
    }
    bool decode(){
        if (!(buffer >> 1u)) buffer = 0x100u + *current++;
        bool bit = buffer & 1u;
        buffer >>= 1u;
        return bit;
    }
    size_t size(){ return (buffer == 1u) ? current - start : current - start + 1; }
private:
    unsigned char const * start = NULL;
    unsigned char const * current = NULL;
    size_t capacity = 0;
    unsigned int buffer = 1u;
};

// size of segment: default 4 MB
const int seg_size = 1;//4 * 1024 * 1024;

// encode the intra level components progressively
/*
@params data: coefficient data
@params n: number of coefficients in current level
@params level_exp: exponent of max level element
*/
template <class T>
vector<unsigned char*> progressive_encoding(T const * data, size_t n, int level_exp){
    vector<unsigned char*> intra_level_components;
    // use 32 bitplanes, where the first one is used for sign
    const int num_bitplanes = 31;
    int num_level_component = n * sizeof(T) / seg_size;
    int consecutive_encoding_bits = 1;
    if(num_level_component >= num_bitplanes) num_level_component = num_bitplanes + 1;
    else{
        // compute consecutive_encoding_bits (CEB)
        num_level_component = 1 << (unsigned int)log2(num_level_component);
        consecutive_encoding_bits = 32 / num_level_component;
    }
    size_t level_component_size = (n * sizeof(T) - 1) / num_level_component + 1;
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << ", consecutive_encoding_bits = " << consecutive_encoding_bits << endl;
    cout << "level_component_size = " << level_component_size << endl;
    vector<BitEncoder> encoders;
    for(int i=0; i<num_level_component; i++){
        unsigned char * buffer = (unsigned char *) malloc(level_component_size);
        intra_level_components.push_back(buffer);
        encoders.push_back(BitEncoder(buffer, level_component_size));
    }
    for(int i=0; i<n; i++){
        T cur_data = ldexp(data[i], num_bitplanes - level_exp);
        long int fix_point = (long int) cur_data;
        // encode each bit of the data for each level component
        // encode sign and the first (CEB - 1) bits
        bool sign = data[i] < 0;
        unsigned int fp = sign ? -fix_point : +fix_point;
        encoders[0].encode(sign);
        unsigned int bit = 31;
        for(int k=1; k<consecutive_encoding_bits; k++){
            encoders[0].encode((fp >> bit) & 1);
            bit --;
        }
        for(int j=1; j<num_level_component; j++){
            // encode CEB bits at a time
            for(int k=0; k<consecutive_encoding_bits; k++){
                encoders[j].encode((fp >> bit) & 1);
                bit --;
            }
        }
    }
    for(int i=0; i<num_level_component; i++){
        // flush current encoding bits
        encoders[i].flush();
    }
    return intra_level_components;
}

template <class T>
T * progressive_decoding(const vector<unsigned char*>& level_components, size_t n, int level_exp, int recompose_level_intra=0){
    T * level_data = (T *) malloc(n * sizeof(T));
    // use 32 bitplanes, where the first one is used for sign
    const int num_bitplanes = 31;
    int num_level_component = n * sizeof(T) / seg_size;
    int consecutive_encoding_bits = 1;
    if(num_level_component >= num_bitplanes) num_level_component = num_bitplanes + 1;
    else{
        // compute consecutive_encoding_bits (CEB)
        num_level_component = 1 << (unsigned int)log2(num_level_component);
        consecutive_encoding_bits = 32 / num_level_component;
    }
    // reconstruct part of the level if requested
    if((recompose_level_intra != 0) && (recompose_level_intra < num_level_component)) num_level_component = recompose_level_intra;
    size_t level_component_size = (n * sizeof(T) - 1) / num_level_component + 1;
    vector<BitDecoder> decoders;
    for(int i=0; i<num_level_component; i++){
        decoders.push_back(BitDecoder(level_components[i], level_component_size));
    }
    T * data_pos = level_data;
    for(int i=0; i<n; i++){
        // decode each bit of the data for each level component
        // decode sign and the first (CEB - 1) bits
        bool sign = decoders[0].decode();
        unsigned int fp = 0;
        unsigned int bit = 31;
        for(int k=1; k<consecutive_encoding_bits; k++){
            unsigned int current_bit = decoders[0].decode();
            fp |= current_bit << bit;
            bit --;
        }
        for(int j=1; j<num_level_component; j++){
            // encode CEB bits at a time
            for(int k=0; k<consecutive_encoding_bits; k++){
                unsigned int current_bit = decoders[j].decode();
                fp |= current_bit << bit;
                bit --;
            }
        }
        signed int fix_point = fp;
        if(sign) fix_point = -fix_point;
        *data_pos = ldexp((float)fix_point, - num_bitplanes + level_exp);
        data_pos ++;
    }
    return level_data;
}

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

// refactor level-centric decomposed data in hierarchical fashion
/*
@params data: decomposed data
@params target_level: decomposed level
@params dims: data dimensions
@params metadata: malloced metadata array
*/
template <class T>
vector<vector<unsigned char*>> level_centric_data_refactor(const T * data, int target_level, const vector<size_t>& dims, unsigned char * metadata){
    int max_level = log2(*min_element(dims.begin(), dims.end()));
    if(target_level > max_level) target_level = max_level;
    auto level_dims = init_levels(dims, target_level);
    size_t * level_elements = reinterpret_cast<size_t*>(metadata);
    T * level_error_bounds = reinterpret_cast<T*>(metadata + (target_level + 1) * sizeof(size_t));
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
    // record all level components
    vector<vector<unsigned char*>> level_components;
    vector<size_t> dims_dummy(dims.size(), 0);
    unsigned char * buffer = (unsigned char *) malloc(level_elements[0] * sizeof(T));
    interleave_level_coefficients(data, dims, level_dims[0], dims_dummy, reinterpret_cast<T*>(buffer), level_error_bounds[0]);
    vector<unsigned char*> first_level;
    first_level.push_back(buffer);
    // first component is the coarsest representation of full precision
    level_components.push_back(first_level);
    for(int i=1; i<=target_level; i++){
        cout << "encoding level " << i << endl;
        unsigned char * buffer = (unsigned char *) malloc(level_elements[i] * sizeof(T));
        // extract components for each level
        interleave_level_coefficients(data, dims, level_dims[i], level_dims[i - 1], reinterpret_cast<T*>(buffer), level_error_bounds[i]);
        if(level_elements[i] * sizeof(T) < seg_size){
            vector<unsigned char*> tiny_level;
            tiny_level.push_back(buffer);
            level_components.push_back(tiny_level);
        }
        else{
            // identify exponent of max element
            int level_exp = 0;
            frexp(level_error_bounds[i], &level_exp);
            // intra-level progressive encoding
            auto intra_level_components = progressive_encoding(reinterpret_cast<T*>(buffer), level_elements[i], level_exp);
            level_components.push_back(intra_level_components);
            // release extracted component
            free(buffer);
        }
    }
    return level_components;
}

// reposition level-centric decomposed data for recomposition
/*
@params data: decomposed data
@params target_level: decomposed level
@params target_recompose_level: recomposed level
@params dims: data dimensions, modified after calling the function
*/
template <class T>
T * level_centric_data_reposition(const vector<vector<unsigned char*>>& level_components, const unsigned char * metadata, int target_level, int target_recompose_level, const vector<int>& intra_recompose_level, vector<size_t>& dims){
    auto level_dims = init_levels(dims, target_level);
    const size_t * level_elements = reinterpret_cast<const size_t*>(metadata);
    const T * level_error_bounds = reinterpret_cast<const T*>(metadata + (target_level + 1) * sizeof(size_t));    
    size_t num_elements = 1;
    for(int j=0; j<dims.size(); j++){
        dims[j] = level_dims[target_recompose_level][j];
        num_elements *= dims[j];
    }
    T * data = (T *) malloc(num_elements * sizeof(T));
    vector<size_t> dims_dummy(dims.size(), 0);
    reposition_level_coefficients(reinterpret_cast<T*>(level_components[0][0]), dims, level_dims[0], dims_dummy, data);
    for(int i=1; i<=target_recompose_level; i++){
        if(level_elements[i] * sizeof(T) < seg_size){
            reposition_level_coefficients(reinterpret_cast<T*>(level_components[i][0]), dims, level_dims[i], level_dims[i - 1], data);
        }
        else{
            // identify exponent of max element
            int level_exp = 0;
            frexp(level_error_bounds[i], &level_exp);
            // intra-level progressive decoding
            auto buffer = progressive_decoding<T>(level_components[i], level_elements[i], level_exp, intra_recompose_level[i]);
            reposition_level_coefficients(buffer, dims, level_dims[i], level_dims[i - 1], data);
            // release reconstructed component
            free(buffer);
        }

    }
    return data;
}

}

#endif