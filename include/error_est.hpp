#ifndef _REFACTOR_ERROR_EST_HPP
#define _REFACTOR_ERROR_EST_HPP

#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <bitset>
#include <queue>
#include "utils.hpp"

namespace REFACTOR{

using namespace std;

#define MAX_ERROR 1
#define SQUARED_ERROR 2
#define PSNR 3

#define MAX_2(a, b) (a > b) ? (a) : (b)
#define LOOKUP_STEPS 1
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
// compute max_e indicator in level
/*
@params data: level data
@params n: number of level data points
@params num_bitplanes: number of encoded bitplanes (include sign)
@params level_exp: aligned level exponent
*/
template <class T>
vector<double> record_level_max_e(const T * data, size_t n, int num_bitplanes, T level_max_val);
template <>
vector<double> record_level_max_e(const float * data, size_t n, int num_bitplanes, float level_max_val){
    int level_exp = 0;
    frexp(level_max_val, &level_exp);
    const int prec = 23;
    int encode_prec = num_bitplanes - 1;
    vector<double> max_e = vector<double>(encode_prec + 1, 0);
    max_e[0] = level_max_val;
    double err = ldexp(1.0, level_exp - 1);
    for(int i=1; i<max_e.size(); i++){
        max_e[i] = err;
        err /= 2;
    }
    return max_e;
}
template <>
vector<double> record_level_max_e(const double * data, size_t n, int num_bitplanes, double level_max_val){
    cout << "Not implemented yet...\nExit -1.\n";
    exit(-1);
    if(num_bitplanes > 52) num_bitplanes = 52;
    vector<double> mse = vector<double>(num_bitplanes, 0);
    // FloatingInt64 fi;
    return mse;
}
// compute mse indicator in level
/*
@params data: level data
@params n: number of level data points
@params num_bitplanes: number of encoded bitplanes (include sign)
@params level_exp: aligned level exponent
*/
template <class T>
vector<double> record_level_mse(const T * data, size_t n, int num_bitplanes, int level_exp);
template <>
vector<double> record_level_mse(const float * data, size_t n, int num_bitplanes, int level_exp){
    const int prec = 23;
    const int encode_prec = num_bitplanes - 1;
    vector<double> mse = vector<double>(num_bitplanes, 0);
    FloatingInt32 fi;
    for(int i=0; i<n; i++){
        if(data[i] == 0) continue;
        int data_exp = 0;
        frexp(data[i], &data_exp);
        auto val = data[i];
        fi.f = val;
        int exp_diff = level_exp - data_exp + prec - encode_prec;
        int index = encode_prec;
        if(exp_diff > 0){
            // zeroing out unrecorded bitplanes
            for(int b=0; b<exp_diff; b++){
                fi.i &= ~(1u << b);            
            }
        }
        else{
            // skip padding 0s (no errors)
            index += exp_diff;
            exp_diff = 0;
        }
        if(index > 0){
            for(int b=exp_diff; b<prec; b++){
                // change b-th bit to 0
                fi.i &= ~(1u << b);
                mse[index] += (data[i] - fi.f)*(data[i] - fi.f);
                index --;
            }
            while(index >= 0){
                mse[index] += data[i] * data[i];
                index --;
            }
        }
    }
    return mse;
}
template <>
vector<double> record_level_mse(const double * data, size_t n, int num_bitplanes, int level_exp){
    cout << "Not implemented yet...\nExit -1.\n";
    exit(-1);
    if(num_bitplanes > 52) num_bitplanes = 52;
    vector<double> mse = vector<double>(num_bitplanes, 0);
    // FloatingInt64 fi;
    return mse;
}

struct Efficiency{
    double efficiency;
    int level;
    Efficiency(double e, int l) : efficiency(e), level(l) {}
};
struct CompareEfficiency { 
    bool operator()(Efficiency const& e1, Efficiency const& e2) { 
        return e1.efficiency < e2.efficiency; 
    } 
}; 

// direct concat the refactored data together
/*
@params level_components: level bitplanes
@params level_sizes: level encoded sizes
@params order: order of level bitplanes
@params total_size: total size of reorganized data
*/
unsigned char * refactored_data_reorganization_in_order(const vector<vector<unsigned char*>>& level_components, const vector<vector<size_t>>& level_sizes, vector<int>& order, size_t& total_size){
    cout << "Reorganize refactored data by level ordering." << endl;
    const int num_levels = level_sizes.size();
    total_size = 0;
    for(int i=0; i<num_levels; i++){
        for(int j=0; j<level_sizes[i].size(); j++){
            total_size += level_sizes[i][j];
        }
    }
    unsigned char * reorganized_data = (unsigned char *) malloc(total_size);
    unsigned char * reorganized_data_pos = reorganized_data;
    for(int i=0; i<num_levels; i++){
        for(int j=0; j<level_sizes[i].size(); j++){
            if(j) order.push_back(i);
            memcpy(reorganized_data_pos, level_components[i][j], level_sizes[i][j]);
            reorganized_data_pos += level_sizes[i][j];
        }
    }
    cout << "recorded data size = " << reorganized_data_pos - reorganized_data << endl;
    for(int i=0; i<order.size(); i++){
        cout << num_levels - 1 - order[i];
    }
    cout << endl;
    return reorganized_data;
}

// organize the refactored data by round-robin
/*
@params level_components: level bitplanes
@params level_sizes: level encoded sizes
@params order: order of level bitplanes
@params total_size: total size of reorganized data
*/
unsigned char * refactored_data_reorganization_round_robin(const vector<vector<unsigned char*>>& level_components, const vector<vector<size_t>>& level_sizes, vector<int>& order, size_t& total_size){
    cout << "Reorganize refactored data by round-robin." << endl;
    const int num_levels = level_sizes.size();
    total_size = 0;
    for(int i=0; i<num_levels; i++){
        for(int j=0; j<level_sizes[i].size(); j++){
            total_size += level_sizes[i][j];
        }
    }
    unsigned char * reorganized_data = (unsigned char *) malloc(total_size);
    unsigned char * reorganized_data_pos = reorganized_data;
    int max_level_size = 0;
    for(int i=0; i<num_levels; i++){
        if(level_sizes[i].size() > max_level_size){
            max_level_size = level_sizes[i].size();
        }
    }
    {
        // sign and first bitplane 
        for(int i=0; i<num_levels; i++){
            order.push_back(i);
            memcpy(reorganized_data_pos, level_components[i][0], level_sizes[i][0]);
            reorganized_data_pos += level_sizes[i][0];
            memcpy(reorganized_data_pos, level_components[i][1], level_sizes[i][1]);
            reorganized_data_pos += level_sizes[i][1];
        }
        for(int j=1; j<max_level_size - 1; j++){
            for(int i=0; i<num_levels; i++){
                if(j >= level_sizes[i].size() - 1) continue;
                order.push_back(i);
                memcpy(reorganized_data_pos, level_components[i][j + 1], level_sizes[i][j + 1]);
                reorganized_data_pos += level_sizes[i][j + 1];
            }
        }
    }
    cout << "recorded data size = " << reorganized_data_pos - reorganized_data << endl;
    for(int i=0; i<order.size(); i++){
        cout << num_levels - 1 - order[i];
    }
    cout << endl;
    return reorganized_data;
}

template <class T>
vector<vector<T>> amortize_error(const vector<vector<T>>& error, int steps){
    vector<vector<T>> amortized = vector<vector<T>>(error.size());
    for(int i=0; i<amortized.size(); i++){
        amortized[i] = vector<T>(error[i].size());
    }
    for(int i=0; i<error.size(); i++){
        for(int j=0; j<error[i].size(); j++){
            double sum = 0;
            int count = 0;
            for(int k=0; k<steps; k++){
                if(j + k == error[i].size()) break;
                sum += error[i][j + k];
                count ++;
            }
            amortized[i][j] = sum / count;
        }
    }
    return amortized;
}

// reorganize refactored data using "uniformed" quantization across levels 
/*
@params N: number of dimensions
@params mode: error estimation mode
@params level_components: level bitplanes
@params level_sizes: level encoded sizes
@params level_errors: level error estimators
@params order: order of level bitplanes
@params total_size: total size of reorganized data
*/
unsigned char * refactored_data_reorganization_uniform_error(int N, int mode, const vector<vector<unsigned char*>>& level_components, const vector<vector<size_t>>& level_sizes, const vector<vector<double>>& level_errors, vector<int>& order, size_t& total_size){
    cout << "Reorganize refactored data by uniform quantization." << endl;
    const int num_levels = level_sizes.size();
    total_size = 0;
    vector<double> factor(num_levels, 1);
    if(mode == SQUARED_ERROR){
        for(int i=0; i<num_levels; i++){
            factor[i] = sqrt(1u << (N * (num_levels - 1 - i)));
            cout << factor[i] << " " << endl;
        }
    }
    cout << "compute total_size" << endl;
    for(int i=0; i<num_levels; i++){
        for(int j=0; j<level_sizes[i].size(); j++){
            total_size += level_sizes[i][j];
        }
    }
    cout << "total_size = " << total_size << endl;
    unsigned char * reorganized_data = (unsigned char *) malloc(total_size);
    unsigned char * reorganized_data_pos = reorganized_data;
    vector<size_t> index(num_levels, 0);
    order.clear();
    int rest_level_count = num_levels;
    double error_level = level_errors[0][0];
    while(rest_level_count){
        for(int i=0; i<num_levels; i++){
            if(index[i] < level_sizes[i].size() - 1){
                while((index[i] < level_sizes[i].size() - 1) && (level_errors[i][index[i]] * factor[i] >= error_level)){
                    order.push_back(i);
                    int level = i;
                    int bitplane_index = index[i];
                    if(bitplane_index == 0){
                        memcpy(reorganized_data_pos, level_components[level][0], level_sizes[level][0]);
                        reorganized_data_pos += level_sizes[level][0];
                        memcpy(reorganized_data_pos, level_components[level][1], level_sizes[level][1]);
                        reorganized_data_pos += level_sizes[level][1];
                    }
                    else{
                        memcpy(reorganized_data_pos, level_components[level][bitplane_index + 1], level_sizes[level][bitplane_index + 1]);
                        reorganized_data_pos += level_sizes[level][bitplane_index + 1];
                    }
                    index[i] ++;
                }
                if(index[i] == level_sizes[i].size() - 1) rest_level_count --;
            }
        }
        error_level /= 2;
    }
    cout << "recorded data size = " << reorganized_data_pos - reorganized_data << endl;
    for(int i=0; i<order.size(); i++){
        cout << num_levels - 1 - order[i];
    }
    cout << endl;
    return reorganized_data;
}

// reorganize refactored data using a greedy estimation
/*
@params N: number of dimensions
@params mode: error estimation mode
@params level_components: level bitplanes
@params level_sizes: level encoded sizes
@params level_errors: level error estimators
@params order: order of level bitplanes
@params total_size: total size of reorganized data
*/
unsigned char * refactored_data_reorganization_greedy_shuffling(int N, int mode, const vector<vector<unsigned char*>>& level_components, const vector<vector<size_t>>& level_sizes, const vector<vector<double>>& level_errors, vector<int>& order, size_t& total_size){
    cout << "Reorganize refactored data by greedy shuffling." << endl;
    const int num_levels = level_sizes.size();
    total_size = 0;
    // init error_gain: reduced error by including current bitplane
    // init sizes: the corresponding sizes of bitplanes
    // NOTE: the sizes of the two vector is 1 less than component sizes
    //         because sign is grouped with the first bitplane
    vector<vector<double>> efficiency;
    for(int i=0; i<num_levels; i++){
        efficiency.push_back(vector<double>(level_sizes[i].size() - 1));
    }
    // compute erorr gain and sizes for bitplanes
    if(mode == MAX_ERROR){
        for(int i=0; i<num_levels; i++){
            for(int j=0; j<efficiency[i].size(); j++){
                auto error_gain = level_errors[i][j] - level_errors[i][j+1];
                auto sizes = (j == 0) ? (level_sizes[i][0] + level_sizes[i][1]) : level_sizes[i][j+1];
                efficiency[i][j] = error_gain / sizes;
                total_size += sizes;
            }
        }
    }
    else if(mode == SQUARED_ERROR){
        vector<double> factor(num_levels, 0);
        for(int i=0; i<num_levels; i++){
            factor[i] = 1u << (N * (num_levels - 1 - i));
        }
        for(int i=0; i<num_levels; i++){
            for(int j=0; j<efficiency[i].size(); j++){
                auto error_gain = (level_errors[i][j] - level_errors[i][j+1]) * factor[i];
                auto sizes = (j == 0) ? (level_sizes[i][0] + level_sizes[i][1]) : level_sizes[i][j+1];
                efficiency[i][j] = error_gain / sizes;
                total_size += sizes;
            }
        }
    }
    else{
        cerr << "Mode " << mode << " not supported in refactored_data_reorganization_shuffled." << endl;
        exit(0);
    }
    // efficiency = amortize_error(efficiency, LOOKUP_STEPS);
    unsigned char * reorganized_data = (unsigned char *) malloc(total_size);
    unsigned char * reorganized_data_pos = reorganized_data;
    vector<size_t> index(num_levels, 0);
    // metric for greedy algorithm: how much error gain per byte
    priority_queue<Efficiency, vector<Efficiency>, CompareEfficiency> efficiency_heap;
    for(int i=0; i<num_levels; i++){
        efficiency_heap.push(Efficiency(efficiency[i][0], i));
    }
    order.clear();
    while(!efficiency_heap.empty()){
        auto eff = efficiency_heap.top();
        efficiency_heap.pop();
        auto level = eff.level;
        auto bitplane_index = index[level];
        // cout << "Encode level " << level << " component " << bitplane_index << ", efficiency = " << eff.efficiency << endl;
        order.push_back(level);
        if(bitplane_index == 0){
            memcpy(reorganized_data_pos, level_components[level][0], level_sizes[level][0]);
            reorganized_data_pos += level_sizes[level][0];
            memcpy(reorganized_data_pos, level_components[level][1], level_sizes[level][1]);
            reorganized_data_pos += level_sizes[level][1];
        }
        else{
            memcpy(reorganized_data_pos, level_components[level][bitplane_index + 1], level_sizes[level][bitplane_index + 1]);
            reorganized_data_pos += level_sizes[level][bitplane_index + 1];
        }
        index[level] ++;
        if(index[level] != level_components[level].size() - 1){
            efficiency_heap.push(Efficiency(efficiency[level][index[level]], level));
        }
    }
    cout << "recorded data size = " << reorganized_data_pos - reorganized_data << endl;
    for(int i=0; i<order.size(); i++){
        cout << num_levels - 1 - order[i];
    }
    cout << endl;
    // exit(0);
    return reorganized_data;
}

inline double estimate_error(double error, int level, int N, int mode){
    if(mode == MAX_ERROR){
        return error * 2.23; // for 3d (1 + 3sqrt(3)/4)
    }
    else if(mode == SQUARED_ERROR){
        return error * (1 << (level * N));
    }
    else{
        cerr << "Mode " << mode << " is not supported! Exit." << endl;
        exit(0);
    }
}
// interpret how many data to read in order to achieve the required tolerance
/*
@params N: dimensions
@params level_sizes: size of level components
@params level_errors: level error estimators
@params order: order of level bitplanes
@params mode: error mode (e.g. max_e = 0, mse = 1)
@params tolerance: required error tolerance
@params index: record number of bitplane extracted from each level
*/
size_t interpret_reading_size(size_t N, const vector<vector<size_t>>& level_sizes, const vector<vector<double>>& level_errors, const vector<int>& order, int mode, double tolerance, vector<int>& index){
    size_t retrieved_size = 0;
    int count = 0;
    double err = 0;
    int num_levels = level_errors.size();
    for(int i=0; i<num_levels; i++){
        err += estimate_error(level_errors[i][0], num_levels - 1 - i, N, mode);
    }
    // for(int i=0; i<level_errors.size(); i++){
    //     for(int j=0; j<level_errors[i].size(); j++){
    //         cout << j << ":" << level_errors[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    while((err > tolerance) && (count < order.size())){
        int level = order[count ++];
        int bitplane_index = index[level];
        if(bitplane_index == 0){
            retrieved_size += level_sizes[level][0] + level_sizes[level][1];
        }
        else{
            retrieved_size += level_sizes[level][bitplane_index + 1];
        }
        err += estimate_error(level_errors[level][bitplane_index + 1], num_levels - 1 - level, N, mode) - estimate_error(level_errors[level][bitplane_index], num_levels - 1 - level, N, mode); 
        index[level] ++;
    }
    cout << "Target error = " << tolerance << endl;
    cout << "Estimated error = " << err << endl;
    return retrieved_size;
}
// read the refactored into vector of level bitplanes
/*
@params refactored_data: refactored data
@params level_sizes: level encoded sizes
@params order: order of level bitplanes
@params retrieved_size: size of retrieved data
*/
vector<vector<const unsigned char*>> read_reorganized_data(const unsigned char * refactored_data, const vector<vector<size_t>>& level_sizes, const vector<int>& order, size_t retrieved_size){
    vector<vector<const unsigned char*>> level_components;
    for(int i=0; i<level_sizes.size(); i++){
        level_components.push_back(vector<const unsigned char*>());
    }
    const unsigned char * refactored_data_pos = refactored_data;
    vector<int> index(level_sizes.size(), 0);
    int count = 0;
    while(refactored_data_pos - refactored_data < retrieved_size){
        int level = order[count ++];
        int bitplane_index = index[level];
        if(bitplane_index == 0){
            level_components[level].push_back(refactored_data_pos);
            refactored_data_pos += level_sizes[level][0];
            level_components[level].push_back(refactored_data_pos);
            refactored_data_pos += level_sizes[level][1];
        }
        else{
            level_components[level].push_back(refactored_data_pos);
            refactored_data_pos += level_sizes[level][bitplane_index + 1];
        }
        index[level] ++;
    }
    cout << "read data size = " << refactored_data_pos - refactored_data << endl; 
    return level_components;
}

}
#endif