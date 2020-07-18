#ifndef _REFACTOR_DATA_ENC_OPT_HPP
#define _REFACTOR_DATA_ENC_OPT_HPP

#include <vector>
#include <cstdlib>

namespace REFACTOR{

using namespace std;

// index for switching to rle
#define HYBRID_ENCODING_SUF_RLE 25
#define EMBEDDED_ENCODING_PRE_RLE 12
#define EMBEDDED_ENCODING_SUF_RLE 25

class EncoderInterface{
public:
    virtual void encode(bool) = 0;
    virtual void flush() = 0;
    virtual size_t size() = 0;
    virtual unsigned char * save(bool use_lossless) = 0;
};

class DecoderInterface{
public:
    virtual ~DecoderInterface() = default;
    virtual bool decode() = 0;
    virtual void load(const unsigned char *) = 0;
};

// modified from ZFP bitstream
class BitEncoder : public EncoderInterface{
protected:
    // bit reversal table used in encoding
    static unsigned char reverse(unsigned char x){
        static const unsigned char lut[] = {
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
    BitEncoder() = default;
    BitEncoder(unsigned char * array) : start(array), current(array){};
    ~BitEncoder() = default;
    inline void encode(bool bit){
        buffer = (buffer << 1u) + bit;
        if(buffer >= 0x100u){
            *current++ = reverse(buffer - 0x100u);
            buffer = 1u;
        }        
    }
    void flush(){
        while (buffer != 1u) encode(false);
    }
    size_t size(){
        if(lossless) return lossless_size;
        return (buffer == 1u) ? current - start : current - start + 1;
    }
    unsigned char * save(bool use_lossless){
        if(use_lossless){
            unsigned char * lossless_compressed = NULL;
            size_t encode_size = (buffer == 1u) ? current - start : current - start + 1;
            lossless_size = MGARD::sz_lossless_compress(ZSTD_COMPRESSOR, 3, start, encode_size, &lossless_compressed);
            lossless = true;
            return lossless_compressed;
        }
        return start;
    }
private:
    unsigned int buffer = 0;
    unsigned char * const start = NULL;
    unsigned char * current = NULL;
    bool lossless = false;
    size_t lossless_size = 0;
};

// A class to read data bit by bit
// Currently does not consider the case when decoding size is larger than capacity
class BitDecoder : public DecoderInterface{
public:
    BitDecoder() = default;
    BitDecoder(bool use_lossless, size_t size) : lossless(use_lossless), lossless_size(size) {};
    ~BitDecoder(){
        if(lossless) free(lossless_decompressed);
    };
    inline bool decode(){
        if (!(buffer >> 1u)) buffer = 0x100u + *current++;
        bool bit = buffer & 1u;
        buffer >>= 1u;
        return bit;
    }
    void load(const unsigned char * encoded_data){
        if(lossless){
            size_t encoded_size = MGARD::sz_lossless_decompress(ZSTD_COMPRESSOR, encoded_data, lossless_size, &lossless_decompressed);
            start = lossless_decompressed;
        }
        else{
            start = encoded_data;
        }
        current = start;
        buffer = 1u;
    }
private:
    unsigned char const * start = NULL;
    unsigned char const * current = NULL;
    unsigned int buffer = 1u;
    bool lossless = false;
    size_t lossless_size = 0;
    unsigned char * lossless_decompressed = NULL;
};
/*******************************************/

// Runlength encoder
#define RLE_CUTOFF_COUNT 256
class RunlengthEncoder : public EncoderInterface{
public:
    RunlengthEncoder(){}
    inline void encode(bool bit){
        if(lastbit == bit){
            count ++;
            if(count == RLE_CUTOFF_COUNT){
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
    unsigned char * save(bool use_lossless){
        // Huffman
        unsigned char * encoded = (unsigned char *) malloc(length.size() * sizeof(int));
        auto encoded_pos = encoded;
        *reinterpret_cast<size_t*>(encoded_pos) = length.size();
        encoded_pos += sizeof(size_t);
        auto encoder = SZ::HuffmanEncoder<int>();
        encoder.preprocess_encode(length, 2*RLE_CUTOFF_COUNT);
        encoder.save(encoded_pos);
        encoder.encode(length, encoded_pos);
        encoder.postprocess_encode();
        encoded_size = encoded_pos - encoded;
        if(use_lossless){
            unsigned char * lossless_compressed = NULL;
            size_t lossless_size = MGARD::sz_lossless_compress(ZSTD_COMPRESSOR, 3, encoded, encoded_size, &lossless_compressed);
            free(encoded);
            encoded_size = lossless_size;
            return lossless_compressed;
        }
        return encoded;
    }
private:
    vector<int> length;
    bool lastbit = false;
    int count = 0;
    int index = 0;
    size_t encoded_size = 0;
};

class RunlengthDecoder : public DecoderInterface{
public:
    RunlengthDecoder() = default;
    RunlengthDecoder(bool use_lossless, size_t size) : lossless(use_lossless), lossless_size(size){};
    ~RunlengthDecoder() = default;
    inline bool decode(){
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
    void load(const unsigned char * encoded){
        const unsigned char * encoded_pos = encoded;
        unsigned char * lossless_compressed = NULL;
        if(lossless){
            size_t encoded_size = MGARD::sz_lossless_decompress(ZSTD_COMPRESSOR, encoded, lossless_size, &lossless_compressed);
            encoded_pos = lossless_compressed;
        }
        size_t n = *reinterpret_cast<const size_t*>(encoded_pos);
        encoded_pos += sizeof(size_t);
        auto encoder = SZ::HuffmanEncoder<int>();
        size_t remaining_length = INT_MAX;
        encoder.load(encoded_pos, remaining_length);
        length = encoder.decode(encoded_pos, n);
        encoder.postprocess_decode();
        if(lossless) free(lossless_compressed);
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
    bool lossless = false;
    size_t lossless_size = 0;
};

// encode bitplanes by byte
/*
@params data: coefficient data
@params n: number of coefficients in current level
@params level_exp: exponent of max level element
@params num_level_component: number of encoded bitplanes
@params byte_encoders: vector of byte-wise encoder
*/
template <class T>
size_t byte_wise_direct_encoding(const T * data, int n, int level_exp, int num_level_component, vector<unsigned char *>& byte_encoders){
    size_t data_index = 0;
    size_t buffer_index = 0;
    for(int i=0; i<n/8; i++){
        unsigned char tmp[32] = {0};
        for(int j=0; j<8; j++){
            T val = data[data_index ++];
            T cur_data = ldexp(val, num_level_component - 1 - level_exp);
            long int fix_point = (long int) cur_data;
            unsigned int sign = val < 0;
            unsigned int fp = sign ? -fix_point : +fix_point;
            tmp[0] += sign << j;
            for(int k=num_level_component - 2; k>=0; k--){
                tmp[num_level_component - 1 - k] += ((fp >> k) & 1) << j;
            }
        }
        for(int k=0; k<num_level_component; k++){
            byte_encoders[k][buffer_index] = tmp[k];
        }
        buffer_index ++;
    }
    {
        int rest = n % 8;
        unsigned char tmp[32] = {0};
        for(int j=0; j<rest; j++){
            T val = data[data_index ++];
            T cur_data = ldexp(val, num_level_component - 1 - level_exp);
            long int fix_point = (long int) cur_data;
            unsigned int sign = val < 0;
            unsigned int fp = sign ? -fix_point : +fix_point;
            tmp[0] += sign << j;
            for(int k=num_level_component - 2; k>=0; k--){
                tmp[num_level_component - 1 - k] += ((fp >> k) & 1) << j;
            }
        }
        for(int k=0; k<num_level_component; k++){
            byte_encoders[k][buffer_index] = tmp[k];
        }
        buffer_index ++;
    }
    return buffer_index;
}

// encode bitplanes by byte (unrolled version)
/*
@params data: coefficient data
@params n: number of coefficients in current level
@params level_exp: exponent of max level element
@params num_level_component: number of encoded bitplanes
@params byte_encoders: vector of byte-wise encoder
*/
template <class T>
size_t byte_wise_direct_encoding_unrolled(const T * data, int n, int level_exp, int num_level_component, vector<unsigned char *>& byte_encoders){
    size_t data_index = 0;
    size_t buffer_index = 0;
    for(int i=0; i<n/8; i++){
        unsigned char tmp = 0;
        T cur_data;
        long int fix_point;
        bool sign;
        cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
        fix_point = (long int) cur_data;
        sign = cur_data < 0;
        tmp += sign;
        unsigned int fp0 = sign ? -fix_point : +fix_point;
        cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
        fix_point = (long int) cur_data;
        sign = cur_data < 0;
        tmp += sign << 1;
        unsigned int fp1 = sign ? -fix_point : +fix_point;
        cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
        fix_point = (long int) cur_data;
        sign = cur_data < 0;
        tmp += sign << 2;
        unsigned int fp2 = sign ? -fix_point : +fix_point;
        cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
        fix_point = (long int) cur_data;
        sign = cur_data < 0;
        tmp += sign << 3;
        unsigned int fp3 = sign ? -fix_point : +fix_point;
        cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
        fix_point = (long int) cur_data;
        sign = cur_data < 0;
        tmp += sign << 4;
        unsigned int fp4 = sign ? -fix_point : +fix_point;
        cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
        fix_point = (long int) cur_data;
        sign = cur_data < 0;
        tmp += sign << 5;
        unsigned int fp5 = sign ? -fix_point : +fix_point;
        cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
        fix_point = (long int) cur_data;
        sign = cur_data < 0;
        tmp += sign << 6;
        unsigned int fp6 = sign ? -fix_point : +fix_point;
        cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
        fix_point = (long int) cur_data;
        sign = cur_data < 0;
        tmp += sign << 7;
        unsigned int fp7 = sign ? -fix_point : +fix_point;
        byte_encoders[0][buffer_index] = tmp;
        for(int k=num_level_component - 2; k>=0; k--){
            byte_encoders[num_level_component - 1 - k][buffer_index] 
                =   ((fp0 >> k) & 1) + (((fp1 >> k) & 1) << 1) 
                    + (((fp2 >> k) & 1) << 2) + (((fp3 >> k) & 1) << 3)
                    + (((fp4 >> k) & 1) << 4) + (((fp5 >> k) & 1) << 5)
                    + (((fp6 >> k) & 1) << 6) + (((fp7 >> k) & 1) << 7);
        }
        buffer_index ++;
    }
    {
        // leftover
        int rest = n % 8;
        unsigned char tmp[32] = {0};
        for(int j=0; j<rest; j++){
            T val = data[data_index ++];
            T cur_data = ldexp(val, num_level_component - 1 - level_exp);
            long int fix_point = (long int) cur_data;
            unsigned int sign = val < 0;
            unsigned int fp = sign ? -fix_point : +fix_point;
            tmp[0] += sign << j;
            for(int k=num_level_component - 2; k>=0; k--){
                tmp[num_level_component - 1 - k] += ((fp >> k) & 1) << j;
            }
        }
        for(int k=0; k<num_level_component; k++){
            byte_encoders[k][buffer_index] = tmp[k];
        }
        buffer_index ++;
    }
    return buffer_index;
}

template <class T>
T * byte_wise_direct_decoding(const vector<const unsigned char*>& level_components, int n, int level_exp, int num_level_component){
    T * level_data = (T *) malloc(n * sizeof(T));
    size_t level_component_size = (n * sizeof(T) - 1) / num_level_component + 1 + 8;
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    cout << "level_component_size = " << level_component_size << endl;
    size_t buffer_index = 0;
    T * data_pos = level_data;
    for(int i=0; i<n/8; i++){
        unsigned int fp0 = 0, fp1 = 0, fp2 = 0, fp3 = 0, fp4 = 0, fp5 = 0, fp6 = 0, fp7 = 0;
        for(int j=1; j<num_level_component; j++){
            unsigned char cur_byte = level_components[j][buffer_index];
            fp0 = (fp0 << 1) + (cur_byte & 1);
            fp1 = (fp1 << 1) + ((cur_byte >> 1) & 1);
            fp2 = (fp2 << 1) + ((cur_byte >> 2) & 1);
            fp3 = (fp3 << 1) + ((cur_byte >> 3) & 1);
            fp4 = (fp4 << 1) + ((cur_byte >> 4) & 1);
            fp5 = (fp5 << 1) + ((cur_byte >> 5) & 1);
            fp6 = (fp6 << 1) + ((cur_byte >> 6) & 1);
            fp7 = (fp7 << 1) + ((cur_byte >> 7) & 1);
            // this is slower
            // fp0 += (cur_byte & 1) << (num_level_component - 1 - j);
            // fp1 += ((cur_byte >> 1) & 1) << (num_level_component - 1 - j);
            // fp2 += ((cur_byte >> 2) & 1) << (num_level_component - 1 - j);
            // fp3 += ((cur_byte >> 3) & 1) << (num_level_component - 1 - j);
            // fp4 += ((cur_byte >> 4) & 1) << (num_level_component - 1 - j);
            // fp5 += ((cur_byte >> 5) & 1) << (num_level_component - 1 - j);
            // fp6 += ((cur_byte >> 6) & 1) << (num_level_component - 1 - j);
            // fp7 += ((cur_byte >> 7) & 1) << (num_level_component - 1 - j);
        }
        unsigned char sign = level_components[0][buffer_index];
        signed int fix_point;
        fix_point = fp0;
        fix_point = (sign & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        fix_point = fp1;
        fix_point = ((sign >> 1) & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        fix_point = fp2;
        fix_point = ((sign >> 2) & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        fix_point = fp3;
        fix_point = ((sign >> 3) & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        fix_point = fp4;
        fix_point = ((sign >> 4) & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        fix_point = fp5;
        fix_point = ((sign >> 5) & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        fix_point = fp6;
        fix_point = ((sign >> 6) & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        fix_point = fp7;
        fix_point = ((sign >> 7) & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        buffer_index ++;
    }
    {
        // leftover
        int rest = n % 8;
        unsigned int fp[8] = {0};
        for(int j=1; j<num_level_component; j++){
            unsigned char cur_byte = level_components[j][buffer_index];
            for(int r=0; r<rest; r++){
                fp[r] = (fp[r] << 1) + ((cur_byte >> r) & 1);
            }
        }
        unsigned char sign = level_components[0][buffer_index];
        for(int r=0; r<rest; r++){
            signed int fix_point = fp[r];
            fix_point = ((sign >> r) & 1) ? -fix_point : fix_point;
            *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        }
        buffer_index ++;
    }
    return level_data;
}

// encode bitplanes by byte (unrolled version), apply sign postpone which records sign after the first 1
/*
@params data: coefficient data
@params n: number of coefficients in current level
@params level_exp: exponent of max level element
@params num_level_component: number of encoded bitplanes
@params byte_encoders: vector of byte-wise encoder
*/
// template <class T>
// size_t byte_wise_direct_encoding_with_sign_postpone(const T * data, int n, int level_exp, int num_level_component, vector<unsigned char *>& byte_encoders){
//     size_t data_index = 0;
//     size_t buffer_index = 0;
//     for(int i=0; i<n/8; i++){
//         unsigned char tmp = 0;
//         T cur_data;
//         long int fix_point;
//         bool sign0, sign1, sign2, sign3, sign4, sign5, sign6, sign7;
//         cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
//         fix_point = (long int) cur_data;
//         sign0 = cur_data < 0;
//         unsigned int fp0 = sign0 ? -fix_point : +fix_point;
//         cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
//         fix_point = (long int) cur_data;
//         sign1 = cur_data < 0;
//         unsigned int fp1 = sign1 ? -fix_point : +fix_point;
//         cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
//         fix_point = (long int) cur_data;
//         sign2 = cur_data < 0;
//         unsigned int fp2 = sign2 ? -fix_point : +fix_point;
//         cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
//         fix_point = (long int) cur_data;
//         sign3 = cur_data < 0;
//         unsigned int fp3 = sign3 ? -fix_point : +fix_point;
//         cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
//         fix_point = (long int) cur_data;
//         sign4 = cur_data < 0;
//         unsigned int fp4 = sign4 ? -fix_point : +fix_point;
//         cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
//         fix_point = (long int) cur_data;
//         sign5 = cur_data < 0;
//         unsigned int fp5 = sign5 ? -fix_point : +fix_point;
//         cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
//         fix_point = (long int) cur_data;
//         sign6 = cur_data < 0;
//         unsigned int fp6 = sign6 ? -fix_point : +fix_point;
//         cur_data = ldexp(data[data_index ++], num_level_component - 1 - level_exp);
//         fix_point = (long int) cur_data;
//         sign7 = cur_data < 0;
//         unsigned int fp7 = sign7 ? -fix_point : +fix_point;
//         byte_encoders[0][buffer_index] = tmp;
//         for(int k=num_level_component - 2; k>=0; k--){
//             byte_encoders[num_level_component - 1 - k][buffer_index] 
//                 =   ((fp0 >> k) & 1) + (((fp1 >> k) & 1) << 1) 
//                     + (((fp2 >> k) & 1) << 2) + (((fp3 >> k) & 1) << 3)
//                     + (((fp4 >> k) & 1) << 4) + (((fp5 >> k) & 1) << 5)
//                     + (((fp6 >> k) & 1) << 6) + (((fp7 >> k) & 1) << 7);
//         }
//         buffer_index ++;
//     }
//     {
//         // leftover
//         int rest = n % 8;
//         unsigned char tmp[32] = {0};
//         for(int j=0; j<rest; j++){
//             T val = data[data_index ++];
//             T cur_data = ldexp(val, num_level_component - 1 - level_exp);
//             long int fix_point = (long int) cur_data;
//             unsigned int sign = val < 0;
//             unsigned int fp = sign ? -fix_point : +fix_point;
//             tmp[0] += sign << j;
//             for(int k=num_level_component - 2; k>=0; k--){
//                 tmp[num_level_component - 1 - k] += ((fp >> k) & 1) << j;
//             }
//         }
//         for(int k=0; k<num_level_component; k++){
//             byte_encoders[k][buffer_index] = tmp[k];
//         }
//         buffer_index ++;
//     }
//     return buffer_index;
// }

template <class T>
T * direct_hybrid_decoding(const vector<const unsigned char*>& level_components, const vector<size_t>& level_sizes, const vector<unsigned char>& use_lossless, size_t n, int level_exp, int num_level_component, const vector<unsigned char>& bitplane_indicator){
    T * level_data = (T *) malloc(n * sizeof(T));
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    vector<DecoderInterface*> decoders;
    for(int i=0; i<num_level_component; i++){
        switch(bitplane_indicator[i]){
            case 0:{
                decoders.push_back(new BitDecoder(use_lossless[i], level_sizes[i]));
                break;
            }
            case 1:{
                decoders.push_back(new RunlengthDecoder(use_lossless[i], level_sizes[i]));
                break;
            }
            default:{
                cerr << "Only direct encoding (indicator = 0) and RLE (indicator = 1) are supported\n";
                exit(0);
            }
        }
        decoders[i]->load(level_components[i]);
    }
    T * data_pos = level_data;
    for(int i=0; i<n; i++){
        // decode each bit of the data for each level component
        bool sign = decoders[0]->decode();
        unsigned int fp = 0;
        for(int j=1; j<num_level_component; j++){
            unsigned int current_bit = decoders[j]->decode();
            fp = (fp << 1) + current_bit;
        }
        long int fix_point = fp;
        if(sign) fix_point = -fix_point;
        *data_pos = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        data_pos ++;
    }
    for(int i=0; i<num_level_component; i++){
        delete decoders[i];
    }
    return level_data;
}

template <class T>
T * byte_wise_hybrid_decoding(const vector<const unsigned char*>& level_components, size_t n, int level_exp, int num_level_component, const vector<unsigned char>& bitplane_indicator){
    T * level_data = (T *) malloc(n * sizeof(T));
    int rle_switch_index = 0;   // index for which the encoder switches from rle to direct
    for(int i=1; i<bitplane_indicator.size(); i++){
        if(bitplane_indicator[i] == 0){
            rle_switch_index = i;
            break;
        }
    }
    if(rle_switch_index > num_level_component) rle_switch_index = num_level_component;
    const int direct_switch_index = (num_level_component > HYBRID_ENCODING_SUF_RLE) ? HYBRID_ENCODING_SUF_RLE : num_level_component;   // index for which the encoder switches from direct to rle
    // decoders for rle in bp 1~switch
    vector<RunlengthDecoder*> pre_decoders;
    for(int i=1; i<rle_switch_index; i++){
        pre_decoders.push_back(new RunlengthDecoder());
        pre_decoders.back()->load(level_components[i]);
    }
    // decoders for rle after bp HYBRID_ENCODING_SUF_RLE
    vector<RunlengthDecoder*> suf_decoders;
    for(int i=HYBRID_ENCODING_SUF_RLE; i<num_level_component; i++){
        suf_decoders.push_back(new RunlengthDecoder());
        suf_decoders.back()->load(level_components[i]);
    }
    size_t buffer_index = 0;
    T * data_pos = level_data;
    for(int i=0; i<n/8; i++){
        unsigned int fp0 = 0, fp1 = 0, fp2 = 0, fp3 = 0, fp4 = 0, fp5 = 0, fp6 = 0, fp7 = 0;
        for(int j=1; j<rle_switch_index; j++){
            fp0 = (fp0 << 1) + pre_decoders[j-1]->decode();
            fp1 = (fp1 << 1) + pre_decoders[j-1]->decode();
            fp2 = (fp2 << 1) + pre_decoders[j-1]->decode();
            fp3 = (fp3 << 1) + pre_decoders[j-1]->decode();
            fp4 = (fp4 << 1) + pre_decoders[j-1]->decode();
            fp5 = (fp5 << 1) + pre_decoders[j-1]->decode();
            fp6 = (fp6 << 1) + pre_decoders[j-1]->decode();
            fp7 = (fp7 << 1) + pre_decoders[j-1]->decode();
        }
        for(int j=rle_switch_index; j<direct_switch_index; j++){
            unsigned char cur_byte = level_components[j][buffer_index];
            fp0 = (fp0 << 1) + (cur_byte & 1);
            fp1 = (fp1 << 1) + ((cur_byte >> 1) & 1);
            fp2 = (fp2 << 1) + ((cur_byte >> 2) & 1);
            fp3 = (fp3 << 1) + ((cur_byte >> 3) & 1);
            fp4 = (fp4 << 1) + ((cur_byte >> 4) & 1);
            fp5 = (fp5 << 1) + ((cur_byte >> 5) & 1);
            fp6 = (fp6 << 1) + ((cur_byte >> 6) & 1);
            fp7 = (fp7 << 1) + ((cur_byte >> 7) & 1);            
        }
        for(int j=HYBRID_ENCODING_SUF_RLE; j<num_level_component; j++){
            fp0 = (fp0 << 1) + suf_decoders[j-HYBRID_ENCODING_SUF_RLE]->decode();
            fp1 = (fp1 << 1) + suf_decoders[j-HYBRID_ENCODING_SUF_RLE]->decode();
            fp2 = (fp2 << 1) + suf_decoders[j-HYBRID_ENCODING_SUF_RLE]->decode();
            fp3 = (fp3 << 1) + suf_decoders[j-HYBRID_ENCODING_SUF_RLE]->decode();
            fp4 = (fp4 << 1) + suf_decoders[j-HYBRID_ENCODING_SUF_RLE]->decode();
            fp5 = (fp5 << 1) + suf_decoders[j-HYBRID_ENCODING_SUF_RLE]->decode();
            fp6 = (fp6 << 1) + suf_decoders[j-HYBRID_ENCODING_SUF_RLE]->decode();
            fp7 = (fp7 << 1) + suf_decoders[j-HYBRID_ENCODING_SUF_RLE]->decode();
        }
        unsigned char sign = level_components[0][buffer_index];
        signed int fix_point;
        fix_point = fp0;
        fix_point = (sign & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        fix_point = fp1;
        fix_point = ((sign >> 1) & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        fix_point = fp2;
        fix_point = ((sign >> 2) & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        fix_point = fp3;
        fix_point = ((sign >> 3) & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        fix_point = fp4;
        fix_point = ((sign >> 4) & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        fix_point = fp5;
        fix_point = ((sign >> 5) & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        fix_point = fp6;
        fix_point = ((sign >> 6) & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        fix_point = fp7;
        fix_point = ((sign >> 7) & 1) ? -fix_point : fix_point;
        *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        buffer_index ++;
    }
    {
        // leftover
        int rest = n % 8;
        unsigned int fp[8] = {0};
        for(int j=1; j<rle_switch_index; j++){
            for(int r=0; r<rest; r++){
                fp[r] = (fp[r] << 1) + pre_decoders[j-1]->decode();
            }
        }        
        for(int j=rle_switch_index; j<direct_switch_index; j++){
            unsigned char cur_byte = level_components[j][buffer_index];
            for(int r=0; r<rest; r++){
                fp[r] = (fp[r] << 1) + ((cur_byte >> r) & 1);
            }
        }
        for(int j=direct_switch_index; j<num_level_component; j++){
            for(int r=0; r<rest; r++){
                fp[r] = (fp[r] << 1) + suf_decoders[j-HYBRID_ENCODING_SUF_RLE]->decode();
            }
        }
        unsigned char sign = level_components[0][buffer_index];
        for(int r=0; r<rest; r++){
            signed int fix_point = fp[r];
            fix_point = ((sign >> r) & 1) ? -fix_point : fix_point;
            *(data_pos ++) = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        }
        buffer_index ++;        
    }
    for(int i=0; i<pre_decoders.size(); i++){
        delete pre_decoders[i];
    }
    for(int i=0; i<suf_decoders.size(); i++){
        delete suf_decoders[i];
    }
    return level_data;
}

}
#endif
