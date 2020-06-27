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
// use ZFP-0.5.5 bitstream
#include "bitstream.h"

namespace REFACTOR{

using namespace std;
using namespace MGARD;

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
    int encoded_bitplanes = 32;

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
    vector<unsigned char *> byte_encoders;
    for(int i=0; i<32; i++){
        unsigned char * buffer = (unsigned char *) malloc(level_component_size);
        intra_level_components.push_back(buffer);
        byte_encoders.push_back(buffer);
    }    
    int index_data = 0;
    int index_buffer = 0;
    for(int i=0; i<n/8; i++){
        unsigned char tmp[32] = {0};
        // unsigned int tmp_fp[8] = {0};
        for(int j=0; j<8; j++){
            T val = data[index_data ++];
            T cur_data = ldexp(val, num_level_component - 1 - level_exp);
            long int fix_point = (long int) cur_data;
            unsigned int sign = val < 0;
            unsigned int fp = sign ? -fix_point : +fix_point;
            tmp[0] += sign << j;
            for(int k=30; k>=0; k--){
                tmp[31 - k] += ((fp >> k) & 1) << j;
            }
        }
        for(int k=0; k<32; k++){
            byte_encoders[k][index_buffer] = tmp[k];
        }
        index_buffer ++;
    }
    {
        int rest = n % 8;
        unsigned char tmp[32] = {0};
        for(int j=0; j<rest; j++){
            T val = data[index_data ++];
            T cur_data = ldexp(val, num_level_component - 1 - level_exp);
            long int fix_point = (long int) cur_data;
            unsigned int sign = val < 0;
            unsigned int fp = sign ? -fix_point : +fix_point;
            tmp[0] += sign << j;
            for(int k=30; k>=0; k--){
                tmp[31 - k] += ((fp >> k) & 1) << j;
            }
        }
        for(int k=0; k<32; k++){
            byte_encoders[k][index_buffer] = tmp[k];
        }
        index_buffer ++;
    }
    for(int k=0; k<32; k++){
        encoded_sizes.push_back(index_buffer);
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
    size_t count = 0;
    for(int i=0; i<num_level_component; i++){
        encoders[i].flush();
        intra_level_components.push_back(encoders[i].save());
        encoded_sizes.push_back(encoders[i].size());
        count += encoders[i].size();
        cout << "Ratio of reading " << i << " bitplanes: " << ((i+1) * (n / 8)) * 1.0 / count << endl;
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
vector<double> record_level_mse(const T * data, size_t n, int num_bitplanes, int level_exp);
template <>
vector<double> record_level_mse(const float * data, size_t n, int num_bitplanes, int level_exp){
    // if(n < 21786949){
    //     return vector<double>();
    // }
    // TODO: bitplanes < 23?
    if(num_bitplanes > 23) num_bitplanes = 23;
    vector<double> mse = vector<double>(num_bitplanes + 1, 0);
    vector<double> max_e = vector<double>(num_bitplanes + 1, 0);
    FloatingInt32 fi;
    for(int i=0; i<n; i++){
        // if(i == 5381319){
        //     cout << setprecision(10) << data[i] << endl;
        // }
        int data_exp = 0;
        frexp(data[i], &data_exp);
        auto val = data[i];
        fi.f = val;
        int exp_diff = level_exp - data_exp;
        // zeroing out unrecorded bitplanes
        for(int b=0; b<exp_diff; b++){
            fi.i &= ~(1u << b);            
        }
        int index = num_bitplanes;
        for(int b=exp_diff; b<num_bitplanes; b++){
            // change b-th bit to 0
            fi.i &= ~(1u << b);
            mse[index] += (data[i] - fi.f)*(data[i] - fi.f);
            // if(i == 5381319){
            //     cout << fi.f << endl;
            // }
            float err = fabs(data[i] - fi.f);
            if(err > max_e[index]) max_e[index] = err;
            index --;
        }
        while(index >= 0){
            mse[index] += data[i] * data[i];
            if(fabs(data[i]) > max_e[index]) max_e[index] = fabs(data[i]);
            index --;
        }
    }
    cout << "\nMAX E in level" << setprecision(4) << endl;
    for(int i=0; i<max_e.size(); i++){
        cout << i << ":" << max_e[i] << " ";
    }
    cout << endl;
    // exit(0);
    return mse;
}
template <>
vector<double> record_level_mse(const double * data, size_t n, int num_bitplanes, int level_exp){
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
    // metadata.mse_estimator = true;
    // record all level components
    vector<vector<unsigned char*>> level_components;
    vector<size_t> dims_dummy(dims.size(), 0);
    for(int i=0; i<=target_level; i++){
        // cout << "encoding level " << i << endl;
        const vector<size_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
        unsigned char * buffer = (unsigned char *) malloc(level_elements[i] * sizeof(T));
        // extract components for each level
        interleave_level_coefficients(data, dims, level_dims[i], prev_dims, reinterpret_cast<T*>(buffer));

        string outfile("decomposed_level_");
        writefile((outfile + to_string(i) + ".dat").c_str(), reinterpret_cast<T*>(buffer), level_elements[i]);

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
            // identify exponent of max element
            int level_exp = 0;
            frexp(level_error_bounds[i], &level_exp);
            cout << "level " << i << " max err = " << level_error_bounds[i] << ", exp = " << level_exp << endl;
            if(metadata.mse_estimator){
                auto level_mse = record_level_mse(reinterpret_cast<T*>(buffer), level_elements[i], metadata.encoded_bitplanes - 1, level_exp);
                metadata.mse.push_back(level_mse);
            }
            // intra-level progressive encoding
            if(metadata.option == ENCODING_DEFAULT){
                // struct timespec start, end;
                // int err = clock_gettime(CLOCK_REALTIME, &start);
                auto intra_level_components = progressive_encoding(reinterpret_cast<T*>(buffer), level_elements[i], level_exp, metadata.encoded_bitplanes, metadata.components_sizes[i]);
                level_components.push_back(intra_level_components);
                // err = clock_gettime(CLOCK_REALTIME, &end);
                // cout << "Byteplane encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
            }
            else if(metadata.option == ENCODING_RLE){
                auto intra_level_components = progressive_encoding_with_rle_compression(reinterpret_cast<T*>(buffer), level_elements[i], level_exp, metadata.encoded_bitplanes, metadata.components_sizes[i]);
                level_components.push_back(intra_level_components);
            }
            else if(metadata.option == ENCODING_LZC){

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
                cout << j << ":" << mse[i][j] << " ";
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
            // cout << "encoded_bitplanes = " << encoded_bitplanes << endl;
            // cout << "MSE: " << i << " " << encoded_bitplanes - 1 << endl; 
            // cout << total_mse << " + " << metadata.mse[i][encoded_bitplanes - 1] << " = ";
            // total_mse += metadata.mse[i][encoded_bitplanes - 1];
            // cout << total_mse << endl;
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
            else if(metadata.option == ENCODING_LZC){

            }
            reposition_level_coefficients(buffer, dims, level_dims[i], prev_dims, data);

            string outfile("reconstructed_level_");
            writefile((outfile + to_string(i) + ".dat").c_str(), reinterpret_cast<T*>(buffer), level_elements[i]);

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