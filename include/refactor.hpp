#ifndef _REFACTOR_HPP
#define _REFACTOR_HPP

#include <vector>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <bitset>
#include "utils.hpp"
// use ZFP-0.5.5 bitstream
#include "bitstream.h"

namespace REFACTOR{

using namespace std;
using namespace MGARD;

// BitEncoder and BitDecoder are modified from ZFP-0.3.1
// Switched to ZFP-0.5.5 for better encoding performance
// Keep original decoder because of performance
/*******************************************/
/*
** Copyright (c) 2014, Lawrence Livermore National Security, LLC.
** Produced at the Lawrence Livermore National Laboratory.
** Written by Peter Lindstrom.
** LLNL-CODE-663824.
** All rights reserved.
*/
#if defined(__GNUC__)
#elif defined(__IBMCPP__)
  #include <builtins.h>
#elif defined(_WIN64)
  #include <intrin.h>
  #ifndef HAVE_C99_MATH
    // for old versions of MSVC that do not have C99 math support
    inline long int lrint(double x) { return  (long int)x; }
    inline long long int llrint(double x) { return (long long int)x; }
  #endif
#else
  #error "compiler not supported"
#endif

inline uint fp_uclz(uint32 x)
{
#if defined(__GNUC__)
  return __builtin_clz(x);
#elif defined(__IBMCPP__)
  return __cntlz4(x);
#elif defined(_WIN64)
  unsigned long n;
  _BitScanReverse(&n, x);
  return 31 - n;
#endif
}

inline uint fp_ufls(uint32 x){
#if defined(__GNUC__) || defined(_WIN64)
  return x ? CHAR_BIT * sizeof(x) - fp_uclz(x) : 0;
#elif defined(__IBMCPP__)
  return CHAR_BIT * sizeof(x) - fp_uclz(x);
#endif
}
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
/*******************************************/
// end ZFP-related

// Runlength encoder
class RunlengthEncoder{
public:
    RunlengthEncoder(){}
    void encode(bool bit){
        if(lastbit == bit){
            count ++;
            if(count == 256){
                length.push_back(count);
                count = 0;
                lastbit = !bit;
            }
        }
        else {
            length.push_back(count);
            count = 1;
            lastbit = bit;
        }
    }
    bool decode(){
        if(count){
            count --;
            return lastbit;
        }
        else{
            count = length[index ++];
            lastbit = !lastbit;
            return decode();
        }
    }
    void flush(){
        if(count != 0){
            length.push_back(count);
            count = 0;
            lastbit = false;
        }
    }
    size_t size(){
        return encoded_size;
    }
    unsigned char * save(){
        // Huffman
        unsigned char * compressed = (unsigned char *) malloc(length.size() * sizeof(int));
        auto compressed_pos = compressed;
        *reinterpret_cast<size_t*>(compressed_pos) = length.size();
        compressed_pos += sizeof(size_t);
        auto encoder = SZ::HuffmanEncoder<int>();
        encoder.preprocess_encode(length, 2*256);
        encoder.save(compressed_pos);
        encoder.encode(length, compressed_pos);
        encoder.postprocess_encode();
        encoded_size = compressed_pos - compressed;
        return compressed;
    }
    void load(const unsigned char * compressed){
        const unsigned char * compressed_pos = compressed;
        size_t n = *reinterpret_cast<const size_t*>(compressed_pos);
        compressed_pos += sizeof(size_t);
        auto encoder = SZ::HuffmanEncoder<int>();
        size_t remaining_length = INT_MAX;
        encoder.load(compressed_pos, remaining_length);
        length = encoder.decode(compressed_pos, n);
        encoder.postprocess_decode();
        count = 0;
        index = 0;
        // toggle lastbit to true such that the first decode would be false since count=0
        lastbit = true;
    }
private:
    vector<int> length;
    bool lastbit = false;
    int count = 0;
    int index = 0;
    int encoded_size = 0;
};

// metadata for refactored data
#define ENCODING_DEFAULT 0
#define ENCODING_RLE 1
#define ENCODING_LZC 2
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
    vector<vector<size_t>> components_sizes; // size of each component in each level
    vector<vector<double>> mse;         // mse[i][j]: mse of level i using the first (j+1) bit-planes
    bool mse_estimator = false;
    int option = ENCODING_DEFAULT;
    int encoded_bitplanes = 24;

    void init_encoded_sizes(){
        components_sizes.clear();
        for(int i=0; i<level_elements.size(); i++){
            components_sizes.push_back(vector<size_t>());
        }
    }
    size_t size(){
        size_t metadata_size = 
                sizeof(int) + sizeof(int) + sizeof(size_t) // option + encoded_bitplanes + number of levels
                + level_elements.size() * sizeof(size_t) 
                + level_error_bounds.size() * sizeof(T)
                + sizeof(char);
        if(mse_estimator) {
            for(int i=0; i<mse.size(); i++){
                metadata_size += sizeof(size_t) + mse[i].size() * sizeof(double);
            }
        }
        return metadata_size;
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
        *buffer_pos = mse_estimator;
        buffer_pos += sizeof(unsigned char);
        if(mse_estimator){
            for(int i=0; i<level_elements.size(); i++){
                *reinterpret_cast<size_t*>(buffer_pos) = mse[i].size();
                buffer_pos += sizeof(size_t);
                memcpy(buffer_pos, mse[i].data(), mse[i].size() * sizeof(double));
                buffer_pos += mse[i].size() * sizeof(double);
            }
        }
        return buffer;
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
        mse_estimator = *buffer_pos;
        buffer_pos += sizeof(unsigned char);
        if(mse_estimator){
            mse.clear();
            for(int i=0; i<level_elements.size(); i++){
                size_t num = *reinterpret_cast<const size_t*>(buffer_pos);
                buffer_pos += sizeof(size_t);
                vector<double> level_mse = vector<double>(reinterpret_cast<const double *>(buffer_pos), reinterpret_cast<const double *>(buffer_pos) + num);
                mse.push_back(level_mse);
                buffer_pos += num * sizeof(double);
            }
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
// size of segment: default 4 MB
const int seg_size = 4 * 1024 * 1024;

// encode the intra level components progressively
/*
@params data: coefficient data
@params n: number of coefficients in current level
@params level_exp: exponent of max level element
*/
template <class T>
vector<unsigned char*> progressive_encoding(T const * data, size_t n, int level_exp, int num_level_component, vector<size_t>& encoded_sizes){
    vector<unsigned char*> intra_level_components;
    size_t level_component_size = (n * sizeof(T) - 1) / num_level_component + 1 + 8;
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    cout << "level_component_size = " << level_component_size << endl;
    // vector<BitEncoder> encoders;
    vector<bitstream*> encoders;
    for(int i=0; i<num_level_component; i++){
        unsigned char * buffer = (unsigned char *) malloc(level_component_size);
        intra_level_components.push_back(buffer);
        // encoders.push_back(BitEncoder(buffer, level_component_size));
        encoders.push_back(stream_open(buffer, level_component_size));
    }
    for(int i=0; i<n; i++){
        T cur_data = ldexp(data[i], num_level_component - 1 - level_exp);
        long int fix_point = (long int) cur_data;
        // encode each bit of the data for each level component
        bool sign = data[i] < 0;
        unsigned int fp = sign ? -fix_point : +fix_point;
        // encoders[0].encode(sign);
        stream_write_bit(encoders[0], sign);
        for(int j=num_level_component - 1; j>0; j--){
            stream_write_bit(encoders[j], fp & 1);
            fp >>= 1;
        }
    }
    for(int i=0; i<num_level_component; i++){
        // flush current encoding bits
        // encoders[i].flush();
        stream_flush(encoders[i]);
        encoded_sizes.push_back(stream_size(encoders[i]));
        stream_close(encoders[i]);
    }
    return intra_level_components;
}

template <class T>
T * progressive_decoding(const vector<unsigned char*>& level_components, size_t n, int level_exp, int num_level_component){
    T * level_data = (T *) malloc(n * sizeof(T));
    size_t level_component_size = (n * sizeof(T) - 1) / num_level_component + 1 + 8;
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    cout << "level_component_size = " << level_component_size << endl;
    vector<BitDecoder> decoders;
    // vector<bitstream*> decoders;
    for(int i=0; i<num_level_component; i++){
        decoders.push_back(BitDecoder(level_components[i], level_component_size));
        // decoders.push_back(stream_open(level_components[i], level_component_size));
    }
    T * data_pos = level_data;
    for(int i=0; i<n; i++){
        // decode each bit of the data for each level component
        bool sign = decoders[0].decode();
        // bool sign = stream_read_bit(decoders[0]);
        unsigned int fp = 0;
        for(int j=1; j<num_level_component; j++){
            unsigned int current_bit = decoders[j].decode();
            // unsigned int current_bit = stream_read_bit(decoders[j]);
            fp = (fp << 1) + current_bit;
        }
        long int fix_point = fp;
        if(sign) fix_point = -fix_point;
        *data_pos = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        data_pos ++;
    }
    // for(int i=0; i<num_level_component; i++){
    //     stream_close(decoders[i]);
    // }
    return level_data;
}

// encode the intra level components progressively, with runlength encoding on each bit-plane
/*
@params data: coefficient data
@params n: number of coefficients in current level
@params level_exp: exponent of max level element
*/
template <class T>
vector<unsigned char*> progressive_encoding_with_rle_compression(T const * data, size_t n, int level_exp, int num_level_component, vector<size_t>& encoded_sizes){
    vector<unsigned char*> intra_level_components;
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    vector<RunlengthEncoder> encoders;
    for(int i=0; i<num_level_component; i++){
        encoders.push_back(RunlengthEncoder());
    }
    for(int i=0; i<n; i++){
        T cur_data = ldexp(data[i], num_level_component - 1 - level_exp);
        long int fix_point = (long int) cur_data;
        // encode each bit of the data for each level component
        bool sign = data[i] < 0;
        unsigned int fp = sign ? -fix_point : +fix_point;
        encoders[0].encode(sign);
        for(int j=num_level_component - 1; j>0; j--){
            encoders[j].encode(fp & 1);
            fp >>= 1;
        }
    }
    for(int i=0; i<num_level_component; i++){
        encoders[i].flush();
        intra_level_components.push_back(encoders[i].save());
        encoded_sizes.push_back(encoders[i].size());
    }
    return intra_level_components;
}

template <class T>
T * progressive_decoding_with_rle_compression(const vector<unsigned char*>& level_components, size_t n, int level_exp, int num_level_component){
    T * level_data = (T *) malloc(n * sizeof(T));
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    vector<RunlengthEncoder> decoders;
    for(int i=0; i<num_level_component; i++){
        decoders.push_back(RunlengthEncoder());
        decoders[i].load(level_components[i]);
    }
    T * data_pos = level_data;
    for(int i=0; i<n; i++){
        // decode each bit of the data for each level component
        bool sign = decoders[0].decode();
        unsigned int fp = 0;
        for(int j=1; j<num_level_component; j++){
            unsigned int current_bit = decoders[j].decode();
            fp = (fp << 1) + current_bit;
        }
        long int fix_point = fp;
        if(sign) fix_point = -fix_point;
        *data_pos = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        data_pos ++;
    }
    return level_data;
}
// encode the intra level components progressively, with lossless compression on leading zeros
/*
@params data: coefficient data
@params n: number of coefficients in current level
@params level_exp: exponent of max level element
*/
template <class T>
vector<unsigned char*> progressive_encoding_with_lzc_compression(T const * data, size_t n, int level_exp, int num_level_component, vector<size_t>& encoded_sizes){
    vector<unsigned char*> intra_level_components;
    size_t level_component_size = (n * sizeof(T) - 1) / num_level_component + 1 + 8;
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    cout << "level_component_size = " << level_component_size << endl;
    const uint prec = 32;
    const uint extra_lzc = prec - num_level_component;
    // TODO: change to unsigned char since max(LZ) < 32?
    vector<int> leading_zeros(n);
    vector<bitstream*> encoders;
    for(int i=0; i<num_level_component; i++){
        unsigned char * buffer = (unsigned char *) malloc(level_component_size);
        encoders.push_back(stream_open(buffer, level_component_size));
    }
    for(int i=0; i<n; i++){
        T cur_data = ldexp(data[i], num_level_component - 1 - level_exp);
        long int fix_point = (long int) cur_data;
        // encode each bit of the data for each level component
        // encode sign
        bool sign = data[i] < 0;
        unsigned int fp = sign ? -fix_point : +fix_point;
        stream_write_bit(encoders[0], sign);
        // record number of leading zeros
        leading_zeros[i] = fp_uclz(fp) - extra_lzc;
        for(int j=num_level_component - 1; j>=leading_zeros[i]; j--){
            stream_write_bit(encoders[j], fp & 1);
            fp >>= 1;
        }
    }
    // record level component size, bitplane size and huffman tree/codes
    unsigned char * bitplane_metadata = (unsigned char *) malloc(n);
    unsigned char * bitplane_metadata_pos = bitplane_metadata;
    auto encoder = SZ::HuffmanEncoder<int>();
    encoder.preprocess_encode(leading_zeros, 2*32);
    encoder.save(bitplane_metadata_pos);
    encoder.encode(leading_zeros, bitplane_metadata_pos);
    encoder.postprocess_encode();
    size_t compressed_length = bitplane_metadata_pos - bitplane_metadata;
    // record lzc component
    encoded_sizes.push_back(compressed_length);
    cout << "Num_elements = " << n << endl;
    cout << "Size after Huffman = " << compressed_length << endl;
    // record huffman size
    // unsigned char * lossless_compressed = NULL;
    // size_t lossless_length = sz_lossless_compress(ZSTD_COMPRESSOR, 3, bitplane_metadata, compressed_length, &lossless_compressed);
    // cout << "Size after lossless = " << lossless_length << endl;
    // size_t count = compressed_length;
    intra_level_components.push_back(bitplane_metadata);
    for(int i=0; i<num_level_component; i++){
        stream_flush(encoders[i]);
        encoded_sizes.push_back(stream_size(encoders[i]));
        intra_level_components.push_back(reinterpret_cast<unsigned char *>(stream_data(encoders[i])));
        stream_close(encoders[i]);
    }
    return intra_level_components;
}

template <class T>
T * progressive_decoding_with_lzc_compression(const vector<unsigned char*>& level_components, size_t n, int level_exp, int num_level_component){
    T * level_data = (T *) malloc(n * sizeof(T));
    size_t level_component_size = (n * sizeof(T) - 1) / num_level_component + 1 + 8;
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    cout << "level_component_size = " << level_component_size << endl;
    // decode leading zeros
    const unsigned char * bitplane_metadata_pos = level_components[0];
    auto encoder = SZ::HuffmanEncoder<int>();
    // remaining_length = n for placeholder (may not be accurate)
    size_t remaining_length = n;
    encoder.load(bitplane_metadata_pos, remaining_length);
    auto leading_zeros = encoder.decode(bitplane_metadata_pos, n);
    encoder.postprocess_decode();
    // decode data
    vector<BitDecoder> decoders;
    for(int i=0; i<num_level_component; i++){
        decoders.push_back(BitDecoder(level_components[i + 1], level_component_size));
    }
    T * data_pos = level_data;
    for(int i=0; i<n; i++){
        // decode each bit of the data for each level component
        // decode sign
        bool sign = decoders[0].decode();
        unsigned int fp = 0;
        for(int j=leading_zeros[i]; j<num_level_component; j++){
            unsigned int current_bit = decoders[j].decode();
            fp = (fp << 1) + current_bit;
        }
        signed int fix_point = fp;
        if(sign) fix_point = -fix_point;
        *data_pos = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
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
void interleave_level_coefficients_3d(const T * data, const vector<size_t>& dims, const vector<size_t>& dims_fine, const vector<size_t>& dims_coasre, T * buffer){
    size_t dim0_offset = dims[1] * dims[2];
    size_t dim1_offset = dims[2];
    size_t count = 0;
    for(int i=0; i<dims_fine[0]; i++){
        for(int j=0; j<dims_fine[1]; j++){
            for(int k=0; k<dims_fine[2]; k++){
                if((i < dims_coasre[0]) && (j < dims_coasre[1]) && (k < dims_coasre[2]))
                    continue;
                buffer[count ++] = data[i*dim0_offset + j*dim1_offset + k];
            }
        }
    }
}
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

// compute maximum value in level
/*
@params data: level data
@params n: number of level data points
*/
template <class T>
T record_level_max_value(const T * data, size_t n){
    T max_val = 0;
    for(int i=0; i<n; i++){
        T val = fabs(data[i]);
        if(val > max_val) max_val = val;
    }
    return max_val;
}

union FloatingInt32{
    float f;
    uint32_t i;
};
union FloatingInt64{
    double f;
    uint64_t i;
};
// compute mse indicator in level
/*
@params data: level data
@params n: number of level data points
@params num_bitplanes: number of encoded bitplanes
*/
template <class T>
vector<double> record_level_mse(const T * data, size_t n, int num_bitplanes);
template <>
vector<double> record_level_mse(const float * data, size_t n, int num_bitplanes){
    // TODO: bitplanes < 23?
    if(num_bitplanes > 23) num_bitplanes = 23;
    vector<double> mse = vector<double>(num_bitplanes, 0);
    FloatingInt32 fi;
    for(int i=0; i<n; i++){
        auto val = data[i];
        fi.f = val;
        for(int b=num_bitplanes - 1; b>=0; b--){
            uint shift = num_bitplanes - 1 - b;
            // change b-th bit to 0
            fi.i &= ~(1u << shift);
            mse[b] += (data[i] - fi.f)*(data[i] - fi.f);
        }
    }
    return mse;
}
template <>
vector<double> record_level_mse(const double * data, size_t n, int num_bitplanes){
    // TODO: bitplanes < 52?
    if(num_bitplanes > 52) num_bitplanes = 52;
    vector<double> mse = vector<double>(num_bitplanes, 0);
    FloatingInt64 fi;
    for(int i=0; i<n; i++){
        auto val = data[i];
        fi.f = val;
        for(int b=num_bitplanes - 1; b>=0; b--){
            uint shift = num_bitplanes - 1 - b;
            // change b-th bit to 0
            fi.i &= ~(1u << shift);
            mse[b] += (data[i] - fi.f)*(data[i] - fi.f);
        }
    }
    return mse;
}
// refactor level-centric decomposed data in hierarchical fashion
/*
@params data: decomposed data
@params target_level: decomposed level
@params dims: data dimensions
@params metadata: malloced metadata array
@params with_compression: whether to use lossless compression on leading zeros
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
    // init metadata sizes recorder
    metadata.init_encoded_sizes();
    // turn on mse
    metadata.mse_estimator = true;
    // record all level components
    vector<vector<unsigned char*>> level_components;
    vector<size_t> dims_dummy(dims.size(), 0);
    for(int i=0; i<=target_level; i++){
        // cout << "encoding level " << i << endl;
        const vector<size_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
        unsigned char * buffer = (unsigned char *) malloc(level_elements[i] * sizeof(T));
        // extract components for each level
        interleave_level_coefficients(data, dims, level_dims[i], prev_dims, reinterpret_cast<T*>(buffer));
        level_error_bounds[i] = record_level_max_value(reinterpret_cast<T*>(buffer), level_elements[i]);
        if(level_elements[i] * sizeof(T) < seg_size){
            if(metadata.mse_estimator){
                auto level_mse = vector<double>(1, 0);
                metadata.mse.push_back(level_mse);
            }
            vector<unsigned char*> tiny_level;
            tiny_level.push_back(buffer);
            level_components.push_back(tiny_level);
            metadata.components_sizes[i].push_back(level_elements[i] * sizeof(T));
        }
        else{
            if(metadata.mse_estimator){
                auto level_mse = record_level_mse(data, level_elements[i], metadata.encoded_bitplanes - 1);
                metadata.mse.push_back(level_mse);
            }
            // identify exponent of max element
            int level_exp = 0;
            frexp(level_error_bounds[i], &level_exp);
            // intra-level progressive encoding
            if(metadata.option == ENCODING_DEFAULT){
                auto intra_level_components = progressive_encoding(reinterpret_cast<T*>(buffer), level_elements[i], level_exp, metadata.encoded_bitplanes, metadata.components_sizes[i]);
                level_components.push_back(intra_level_components);
            }
            else if(metadata.option == ENCODING_RLE){
                auto intra_level_components = progressive_encoding_with_rle_compression(reinterpret_cast<T*>(buffer), level_elements[i], level_exp, metadata.encoded_bitplanes, metadata.components_sizes[i]);
                level_components.push_back(intra_level_components);
            }
            else if(metadata.option == ENCODING_LZC){
                auto intra_level_components = progressive_encoding_with_lzc_compression(reinterpret_cast<T*>(buffer), level_elements[i], level_exp, metadata.encoded_bitplanes, metadata.components_sizes[i]);
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
@params target_level: decomposed level
@params target_recompose_level: recomposed level
@params dims: data dimensions, modified after calling the function
@params with_compression: whether to use lossless compression on leading zeros
*/
template <class T>
T * level_centric_data_reposition(const vector<vector<unsigned char*>>& level_components, const Metadata<T>& metadata, int target_level, int target_recompose_level, const vector<int>& intra_recompose_level, vector<size_t>& dims){
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
    if(metadata.mse_estimator){
        auto& mse = metadata.mse;
        for(int i=0; i<mse.size(); i++){
            for(int j=0; j<mse[i].size(); j++){
                cout << mse[i][j] << " ";
            }
            cout <<endl;
        }
        cout << endl;
    }
    double total_mse = 0;
    // reposition_level_coefficients(reinterpret_cast<T*>(level_components[0][0]), dims, level_dims[0], dims_dummy, data);
    for(int i=0; i<=target_recompose_level; i++){
        const vector<size_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
        if(level_elements[i] * sizeof(T) < seg_size){
            reposition_level_coefficients(reinterpret_cast<T*>(level_components[i][0]), dims, level_dims[i], prev_dims, data);
        }
        else{
            cout << "decoding level " << i << ", size of components = " << level_components[i].size() << endl;
            // identify exponent of max element
            int level_exp = 0;
            frexp(level_error_bounds[i], &level_exp);
            T * buffer = NULL;
            int encoded_bitplanes = intra_recompose_level[i] ? intra_recompose_level[i] : metadata.encoded_bitplanes;
            if(encoded_bitplanes > metadata.encoded_bitplanes) encoded_bitplanes = metadata.encoded_bitplanes;
            total_mse += metadata.mse[i][encoded_bitplanes];
            // intra-level progressive decoding
            if(metadata.option == ENCODING_DEFAULT){
                buffer = progressive_decoding<T>(level_components[i], level_elements[i], level_exp, encoded_bitplanes);                
            }
            else if(metadata.option == ENCODING_RLE){
                buffer = progressive_decoding_with_rle_compression<T>(level_components[i], level_elements[i], level_exp, encoded_bitplanes);
            }
            else if(metadata.option == ENCODING_LZC){
                buffer = progressive_decoding_with_lzc_compression<T>(level_components[i], level_elements[i], level_exp, encoded_bitplanes);
            }
            reposition_level_coefficients(buffer, dims, level_dims[i], prev_dims, data);
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