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
#include <climits>

namespace REFACTOR{

using namespace std;

#define ENCODING_DEFAULT 0
#define ENCODING_DEFAULT_SIGN_POSTPONE 1
#define ENCODING_RLE 2
#define ENCODING_HYBRID 3
#define ENCODING_HYBRID_SIGN_POSTPONE 4

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

template <class T>
vector<unsigned char*> progressive_encoding_with_sign_postpone(T const * data, size_t n, int level_exp, int num_level_component, vector<size_t>& encoded_sizes){
    vector<unsigned char*> intra_level_components;
    size_t level_component_size = (n * sizeof(T) - 1) / num_level_component + 1 + 8;
    cout << "Using direct encoding with sign postpone" << endl;
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    // push back sign bit-plane, nothing in this case
    intra_level_components.push_back(NULL);
    for(int i=1; i<num_level_component; i++){
        unsigned char * buffer = (unsigned char *) malloc(2 * level_component_size);
        intra_level_components.push_back(buffer);
    }
    byte_wise_direct_encoding_with_sign_postpone(data, n, level_exp, num_level_component, intra_level_components, encoded_sizes);
    return intra_level_components;    
}

template <class T>
T * progressive_decoding_with_sign_postpone(const vector<const unsigned char*>& level_components, size_t n, int level_exp, int num_level_component){
    T * level_data = (T *) malloc(n * sizeof(T));
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    vector<BitDecoder*> decoders;
    for(int i=1; i<num_level_component; i++){
        decoders.push_back(new BitDecoder());
        decoders.back()->load(level_components[i]);
    }
    T * data_pos = level_data;
    for(int i=0; i<n; i++){
        // decode each bit of the data for each level component
        bool sign = false;
        unsigned int fp = 0;
        bool first_bit = true;
        for(int j=1; j<num_level_component; j++){
            bool current_bit = decoders[j-1]->decode();
            if(current_bit && first_bit){
                // decode sign
                sign = decoders[j-1]->decode();
                first_bit = false;
            }
            fp = (fp << 1) + current_bit;
        }
        long int fix_point = fp;
        if(sign) fix_point = -fix_point;
        *data_pos = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        data_pos ++;
    }
    for(int i=0; i<decoders.size(); i++){
        delete decoders[i];
    }
    return level_data;
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
        if((k >= HYBRID_ENCODING_SUF_RLE) || use_rle){
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

template <class T>
vector<unsigned char*> progressive_hybrid_embedded_encoding(T const * data, size_t n, int level_exp, int num_level_component, vector<size_t>& encoded_sizes, vector<unsigned char>& bitplane_indicator){
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
    vector<RunlengthEncoder*> pre_encoders;
    for(int i=1; i<EMBEDDED_ENCODING_PRE_RLE; i++){
        pre_encoders.push_back(new RunlengthEncoder());
    }
    // skip sign because sign is recorded with the first value data
    for(int i=0; i<buffer_index; i++){
        unsigned char sign_datum = byte_encoders[0][i];
        bool sign_recorded[8] = {false};
        // iterate bit-planes
        for(int k=1; k<EMBEDDED_ENCODING_PRE_RLE; k++){
            unsigned char datum = byte_encoders[k][i];
            // iterate bit of extracted char
            for(int j=0; j<8; j++){
                if(datum & 1){
                    pre_encoders[k - 1]->encode(true);
                    if(!sign_recorded[j]){
                        // record sign
                        pre_encoders[k - 1]->encode((sign_datum >> j) & 1);
                        sign_recorded[j] = true;
                    }
                }
                else{
                    pre_encoders[k - 1]->encode(false);
                }
                datum >>= 1;
                // if no '1' appears
                // encode to the last RLE bit-plane        
                if((k == EMBEDDED_ENCODING_PRE_RLE - 1) && (!sign_recorded[j])){
                    pre_encoders.back()->encode((sign_datum >> j) & 1);
                }
            }
        }
    }
    // sign is distributed to other bit-planes
    free(intra_level_components[0]);
    intra_level_components[0] = NULL;
    encoded_sizes[0] = 0;
    bitplane_indicator[0] = 1;
    for(int k=1; k<EMBEDDED_ENCODING_PRE_RLE; k++){
        free(intra_level_components[k]);
        pre_encoders[k - 1]->flush();
        intra_level_components[k] = pre_encoders[k - 1]->save();
        encoded_sizes[k] = pre_encoders[k - 1]->size();
        bitplane_indicator[k] = 1;
        delete pre_encoders[k - 1];
    }
    for(int k=EMBEDDED_ENCODING_SUF_RLE; k<num_level_component; k++){
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
    }
    return intra_level_components;
}

template <class T>
T * progressive_hybrid_embedded_decoding(const vector<const unsigned char*>& level_components, size_t n, int level_exp, int num_level_component, const vector<unsigned char>& bitplane_indicator){
    T * level_data = (T *) malloc(n * sizeof(T));
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    vector<DecoderInterface*> decoders;
    const int rle_switch_index_1 = (EMBEDDED_ENCODING_PRE_RLE > num_level_component) ? num_level_component : EMBEDDED_ENCODING_PRE_RLE;
    const int rle_switch_index_2 = (EMBEDDED_ENCODING_SUF_RLE > num_level_component) ? num_level_component : EMBEDDED_ENCODING_SUF_RLE;
    for(int i=1; i<rle_switch_index_1; i++){
        decoders.push_back(new RunlengthDecoder());
        decoders.back()->load(level_components[i]);
    }
    for(int i=rle_switch_index_1; i<rle_switch_index_2; i++){
        decoders.push_back(new BitDecoder());
        decoders.back()->load(level_components[i]);        
    }
    for(int i=rle_switch_index_2; i<num_level_component; i++){
        decoders.push_back(new RunlengthDecoder());
        decoders.back()->load(level_components[i]);
    }
    T * data_pos = level_data;
    for(int i=0; i<n; i++){
        // decode each bit of the data for each level component
        bool sign = false;
        unsigned int fp = 0;
        bool first_bit = true;
        // decode the first rle_switch_index_1 bits to get sign
        for(int j=1; j<rle_switch_index_1; j++){
            unsigned int current_bit = decoders[j-1]->decode();
            if(current_bit && first_bit){
                // decode sign
                sign = decoders[j-1]->decode();
                first_bit = false;
            }
            fp = (fp << 1) + current_bit;
        }
        // if no '1' appears
        if((EMBEDDED_ENCODING_PRE_RLE <= num_level_component) && first_bit){
            sign = decoders[rle_switch_index_1 - 1 - 1]->decode();
        }
        for(int j=rle_switch_index_1; j<num_level_component; j++){
            unsigned int current_bit = decoders[j-1]->decode();
            fp = (fp << 1) + current_bit;            
        }
        long int fix_point = fp;
        if(sign) fix_point = -fix_point;
        *data_pos = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        data_pos ++;
    }
    for(int i=0; i<decoders.size(); i++){
        delete decoders[i];
    }
    return level_data;
}

}
#endif
