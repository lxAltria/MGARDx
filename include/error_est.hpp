#ifndef _REFACTOR_ERROR_EST_HPP
#define _REFACTOR_ERROR_EST_HPP

#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <bitset>

namespace REFACTOR{

using namespace std;

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
    // TODO: bitplanes < 23?
    if(num_bitplanes > 23) num_bitplanes = 23;
    vector<double> mse = vector<double>(num_bitplanes + 1, 0);
    // vector<double> max_e = vector<double>(num_bitplanes + 1, 0);
    FloatingInt32 fi;
    for(int i=0; i<n; i++){
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
            // float err = fabs(data[i] - fi.f);
            // if(err > max_e[index]) max_e[index] = err;
            index --;
        }
        while(index >= 0){
            mse[index] += data[i] * data[i];
            // if(fabs(data[i]) > max_e[index]) max_e[index] = fabs(data[i]);
            index --;
        }
    }
    // cout << "\nMAX E in level" << setprecision(4) << endl;
    // for(int i=0; i<max_e.size(); i++){
    //     cout << i << ":" << max_e[i] << " ";
    // }
    // cout << endl;
    return mse;
}
template <>
vector<double> record_level_mse(const double * data, size_t n, int num_bitplanes, int level_exp){
    cout << "Not implemented yet...\nExit -1.\n";
    exit(-1);
    // TODO: bitplanes < 52?
    if(num_bitplanes > 52) num_bitplanes = 52;
    vector<double> mse = vector<double>(num_bitplanes, 0);
    // FloatingInt64 fi;
    return mse;
}

}
#endif