#ifndef _REFACTOR_ERROR_EST_HPP
#define _REFACTOR_ERROR_EST_HPP

#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <bitset>
#include <queue>

namespace REFACTOR{

using namespace std;

#define MAX_2(a, b) (a > b) ? (a) : (b)
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
vector<double> record_level_max_e(const T * data, size_t n, int num_bitplanes, int level_exp);
template <>
vector<double> record_level_max_e(const float * data, size_t n, int num_bitplanes, int level_exp){
    const int prec = 23;
    int encode_prec = num_bitplanes - 1;
    // if(encode_prec > prec) encode_prec = prec;
    vector<double> max_e = vector<double>(encode_prec + 1, 0);
    FloatingInt32 fi;
    for(int i=0; i<n; i++){
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
        for(int b=exp_diff; b<prec; b++){
            // change b-th bit to 0
            fi.i &= ~(1u << b);
            max_e[index] = MAX_2(max_e[index], fabs(data[i] - fi.f));
            index --;
        }
        while(index >= 0){
            max_e[index] = MAX_2(max_e[index], fabs(data[i]));
            index --;
        }
    }
    return max_e;
}
template <>
vector<double> record_level_max_e(const double * data, size_t n, int num_bitplanes, int level_exp){
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
// reorganize refactored data with respect to error estimators
// use a greedy algorithm to pick up the most efficient bitplane
// return reorganized data
/*
@params level_components: level bitplanes
@params level_sizes: level encoded sizes
@params level_errors: level error estimators
@params order: order of level bitplanes
@params total_size: total size of reorganized data
*/
unsigned char * refactored_data_reorganization(const vector<vector<unsigned char*>>& level_components, const vector<vector<size_t>>& level_sizes, const vector<vector<double>>& level_errors, vector<int>& order, size_t& total_size){
    const int num_levels = level_components.size();
    total_size = 0;
    // init error_gain: reduced error by including current bitplane
    // init sizes: the corresponding sizes of bitplanes
    // NOTE: the sizes of the two vector is 1 less than component sizes
    //         because sign is grouped with the first bitplane
    vector<vector<double>> error_gain;
    vector<vector<double>> sizes;
    for(int i=0; i<num_levels; i++){
        error_gain.push_back(vector<double>(level_components[i].size() - 1));
        sizes.push_back(vector<double>(level_components[i].size() - 1));
    }
    // compute erorr gain and sizes for bitplanes
    for(int i=0; i<num_levels; i++){
        for(int j=0; j<error_gain[i].size(); j++){
            error_gain[i][j] = level_errors[i][j] - level_errors[i][j+1];
            sizes[i][j] = (j == 0) ? (level_sizes[i][0] + level_sizes[i][1]) : level_sizes[i][j+1];
            total_size += sizes[i][j];
        }
    }
    cout << "total_size = " << total_size << endl;
    unsigned char * reorganized_data = (unsigned char *) malloc(total_size);
    unsigned char * reorganized_data_pos = reorganized_data;
    vector<size_t> index(num_levels, 0);
    // metric for greedy algorithm: how much error gain per byte
    priority_queue<Efficiency, vector<Efficiency>, CompareEfficiency> efficiency_heap;
    for(int i=0; i<num_levels; i++){
        efficiency_heap.push(Efficiency(error_gain[i][0] / sizes[i][0], i));
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
            efficiency_heap.push(Efficiency(error_gain[level][bitplane_index] / sizes[level][bitplane_index], level));
        }
    }
    cout << "recorded data size = " << reorganized_data_pos - reorganized_data << endl;
    return reorganized_data;
}

// interpret how many data to read in order to achieve the required tolerance
/*
@params level_sizes: size of level components
@params level_errors: level error estimators
@params order: order of level bitplanes
@params mode: error mode (e.g. max_e = 0, mse = 1)
@params tolerance: required error tolerance
@params index: record number of bitplane extracted from each level
*/
size_t interpret_reading_size(const vector<vector<size_t>>& level_sizes, const vector<vector<double>>& level_errors, const vector<int>& order, int mode, double tolerance, vector<int>& index){
    size_t retrieved_size = 0;
    int count = 0;
    double err = 0;
    for(int i=0; i<level_errors.size(); i++){
        err += level_errors[i][0];
    }
    for(int i=0; i<level_errors.size(); i++){
        for(int j=0; j<level_errors[i].size(); j++){
            cout << j << ":" << level_errors[i][j] << " ";
        }
        cout << endl;
    }
    double err_est_constant = 1.65;
    while((err_est_constant * err > tolerance) && (count < order.size())){
        int level = order[count ++];
        int bitplane_index = index[level];
        if(bitplane_index == 0){
            retrieved_size += level_sizes[level][0] + level_sizes[level][1];
        }
        else{
            retrieved_size += level_sizes[level][bitplane_index + 1];
        }
        err += level_errors[level][bitplane_index + 1] - level_errors[level][bitplane_index]; 
        index[level] ++;
        cout << err << endl;
    }
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