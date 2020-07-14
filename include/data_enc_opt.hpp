#ifndef _REFACTOR_DATA_ENC_OPT_HPP
#define _REFACTOR_DATA_ENC_OPT_HPP

#include <vector>
#include <cstdlib>

namespace REFACTOR{

using namespace std;
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
            // tmp = 0;
            // tmp += (fp0 >> k) & 1;
            // tmp += ((fp1 >> k) & 1) << 1;
            // tmp += ((fp2 >> k) & 1) << 2;
            // tmp += ((fp3 >> k) & 1) << 3;
            // tmp += ((fp4 >> k) & 1) << 4;
            // tmp += ((fp5 >> k) & 1) << 5;
            // tmp += ((fp6 >> k) & 1) << 6;
            // tmp += ((fp7 >> k) & 1) << 7;
            // byte_encoders[num_level_component - 1 - k][buffer_index] = tmp;
        }
        // unsigned char tmp_[32] = {0};
        // for(int k=30; k>=0; k--){
        //     tmp_[31 - k] += (fp0 >> k) & 1;
        // }
        // for(int k=30; k>=0; k--){
        //     tmp_[31 - k] += ((fp1 >> k) & 1) << 1;
        // }
        // for(int k=30; k>=0; k--){
        //     tmp_[31 - k] += ((fp2 >> k) & 1) << 2;
        // }
        // for(int k=30; k>=0; k--){
        //     tmp_[31 - k] += ((fp3 >> k) & 1) << 3;
        // }
        // for(int k=30; k>=0; k--){
        //     tmp_[31 - k] += ((fp4 >> k) & 1) << 4;
        // }
        // for(int k=30; k>=0; k--){
        //     tmp_[31 - k] += ((fp5 >> k) & 1) << 5;
        // }
        // for(int k=30; k>=0; k--){
        //     tmp_[31 - k] += ((fp6 >> k) & 1) << 6;
        // }
        // for(int k=30; k>=0; k--){
        //     tmp_[31 - k] += ((fp7 >> k) & 1) << 7;
        // }
        // for(int k=1; k<32; k++){
        //     byte_encoders[k][buffer_index] = tmp_[k];
        // }
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

}
#endif
