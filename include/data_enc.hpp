#ifndef _REFACTOR_DATA_ENC_HPP
#define _REFACTOR_DATA_ENC_HPP

#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <bitset>
#include <iomanip>
#include "data_enc_opt.hpp"

namespace REFACTOR{

using namespace std;

#define ENCODING_DEFAULT 0
#define ENCODING_RLE 1
#define ENCODING_HYBRID 2

// encode the intra level components progressively
/*
@params data: coefficient data
@params n: number of coefficients in current level
@params level_exp: exponent of max level element
@params num_level_component: number of encoded bitplanes
@params encoded_sizes: size of encoded data
*/
template <class T>
vector<unsigned char*> progressive_encoding(T const * data, size_t n, int level_exp, int num_level_component, vector<size_t>& encoded_sizes){
    vector<unsigned char*> intra_level_components;
    size_t level_component_size = (n * sizeof(T) - 1) / num_level_component + 1 + 8;
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    cout << "level_component_size = " << level_component_size << endl;
    vector<unsigned char *> byte_encoders;
    for(int i=0; i<num_level_component; i++){
        unsigned char * buffer = (unsigned char *) malloc(level_component_size);
        intra_level_components.push_back(buffer);
        byte_encoders.push_back(buffer);
    }    
    size_t buffer_index = byte_wise_direct_encoding_unrolled(data, n, level_exp, num_level_component, byte_encoders);
    for(int k=0; k<num_level_component; k++){
        encoded_sizes.push_back(buffer_index);
    }
    return intra_level_components;
}


template <class T>
T * progressive_decoding(const vector<const unsigned char*>& level_components, size_t n, int level_exp, int num_level_component){
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    return byte_wise_direct_decoding<T>(level_components, n, level_exp, num_level_component);
}

// encode the intra level components progressively, with runlength encoding on each bit-plane
/*
@params data: coefficient data
@params n: number of coefficients in current level
@params level_exp: exponent of max level element
@params num_level_component: number of encoded bitplanes
@params encoded_sizes: size of encoded data
*/
template <class T>
vector<unsigned char*> progressive_encoding_with_rle_compression(T const * data, size_t n, int level_exp, int num_level_component, vector<size_t>& encoded_sizes){
    vector<unsigned char*> intra_level_components;
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    struct timespec start, end;
    int err = clock_gettime(CLOCK_REALTIME, &start);        
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
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "RLE encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    size_t count = 0;
    for(int i=0; i<num_level_component; i++){
        encoders[i].flush();
        err = clock_gettime(CLOCK_REALTIME, &start);        
        intra_level_components.push_back(encoders[i].save());
        err = clock_gettime(CLOCK_REALTIME, &end);
        encoded_sizes.push_back(encoders[i].size());
        count += encoders[i].size();
        // cout << "The " << i << "-th bitplanes encoding time = " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << " s, ratio = " << (n / 8) * 1.0 /encoders[i].size()<< " , progressive ratio = " << ((i+1) * (n / 8)) * 1.0 / count << endl;
    }
    return intra_level_components;
}

template <class T>
T * progressive_decoding_with_rle_compression(const vector<const unsigned char*>& level_components, size_t n, int level_exp, int num_level_component){
    T * level_data = (T *) malloc(n * sizeof(T));
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    vector<RunlengthDecoder> decoders;
    for(int i=0; i<num_level_component; i++){
        decoders.push_back(RunlengthDecoder());
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

// encode the intra level components progressively using alternatives of direct encoding and rle encoding
/*
@params data: coefficient data
@params n: number of coefficients in current level
@params level_exp: exponent of max level element
@params num_level_component: number of encoded bitplanes
@params encoded_sizes: size of encoded data
@params bitplane_indicator: indicator for encoder selection
*/
template <class T>
vector<unsigned char*> progressive_hybrid_encoding(T const * data, size_t n, int level_exp, int num_level_component, vector<size_t>& encoded_sizes, vector<unsigned char>& bitplane_indicator){
    vector<unsigned char*> intra_level_components;
    size_t level_component_size = (n * sizeof(T) - 1) / num_level_component + 1 + 8;
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    cout << "level_component_size = " << level_component_size << endl;
    vector<unsigned char *> byte_encoders;
    for(int i=0; i<num_level_component; i++){
        unsigned char * buffer = (unsigned char *) malloc(level_component_size);
        intra_level_components.push_back(buffer);
        byte_encoders.push_back(buffer);
        // set indicator
        bitplane_indicator.push_back(0);
    }    
    size_t buffer_index = byte_wise_direct_encoding_unrolled(data, n, level_exp, num_level_component, byte_encoders);
    for(int k=0; k<num_level_component; k++){
        encoded_sizes.push_back(buffer_index);
    }
    struct timespec start, end;
    bool use_rle = true;
    for(int k=1; k<num_level_component; k++){
        if((k >= RLE_2_INDEX_F32) || use_rle){
            // int err = clock_gettime(CLOCK_REALTIME, &start);        
            RunlengthEncoder rle;
            for(int i=0; i<buffer_index; i++){
                unsigned char datum = byte_encoders[k][i];
                for(int j=0; j<8; j++){
                    rle.encode(datum & 1);
                    datum >>= 1;            
                }
            }
            rle.flush();
            free(intra_level_components[k]);
            // change content of level components, encoded size and indicator
            intra_level_components[k] = rle.save();
            encoded_sizes[k] = rle.size();
            bitplane_indicator[k] = 1;
            if(rle.size() * 1.5 > level_component_size) use_rle = false;
            // cout << "RLE index = " << k << endl;  
            // err = clock_gettime(CLOCK_REALTIME, &end);
            // cout << "bitplane " << k << " runlength encoding time = " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;
            // cout << "s, encoded size = " << rle.size() << endl;
        }
    }
    return intra_level_components;
}

template <class T>
T * progressive_hybrid_decoding(const vector<const unsigned char*>& level_components, size_t n, int level_exp, int num_level_component, const vector<unsigned char>& bitplane_indicator){
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    return byte_wise_hybrid_decoding<T>(level_components, n, level_exp, num_level_component, bitplane_indicator);
}

}
#endif