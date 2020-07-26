#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <sys/stat.h>
#include "refactor.hpp"
#include "error_est.hpp"

using namespace std;

template <class T>
T * test_reposition(int retrieve_mode, double tolerance, vector<size_t>& recompose_dims, size_t& recompose_times){
    REFACTOR::Metadata<T> metadata;
    metadata.from_file(string("refactor_data/metadata").c_str());
    if(retrieve_mode == PSNR){    
        // change tolerance from psnr to se
        double psnr = tolerance;
        double value_range = metadata.max_val - metadata.min_val;
        size_t num = 1;
        for(const auto& d:metadata.dims){
            num *= d;
        }
        tolerance = (value_range / pow(10, psnr/20))*(value_range / pow(10, psnr/20)) * num;
        metadata.set_mode(SQUARED_ERROR);
    }
    auto data = multigrid_data_recompose(string("refactor_data/refactored.dat"), metadata, tolerance, recompose_dims, recompose_times);
    cout << "Recomposed dims: ";
    size_t num_elements = 1;
    for(int i=0; i<recompose_dims.size(); i++){
        cout << recompose_dims[i] << " ";
        num_elements *= recompose_dims[i];
    }
    cout << endl;
    MGARD::writefile(string("mgard.recomposed").c_str(), data, num_elements);
    return data;
}

// place level-centric data to original data dimensions
template <class T>
void put_back_data_1d(const T * data, const vector<size_t>& dims, T * full_data, const vector<size_t>& full_dims){
    memcpy(full_data, data, dims[0] * sizeof(T));
}
template <class T>
void put_back_data_2d(const T * data, const vector<size_t>& dims, T * full_data, const vector<size_t>& full_dims){
    const T * data_pos = data;
    T * full_data_pos = full_data;
    for(int i=0; i<dims[0]; i++){
        memcpy(full_data_pos, data_pos, dims[1] * sizeof(T));
        data_pos += dims[1];
        full_data_pos += full_dims[1];
    }
}
template <class T>
void put_back_data_3d(const T * data, const vector<size_t>& dims, T * full_data, const vector<size_t>& full_dims){
    const T * data_pos = data;
    T * full_data_pos = full_data;
    for(int i=0; i<dims[0]; i++){
        const T * row_data_pos = data_pos;
        T * row_full_data_pos = full_data_pos;
        for(int j=0; j<dims[1]; j++){
            memcpy(row_full_data_pos, row_data_pos, dims[2] * sizeof(T));
            row_data_pos += dims[2];
            row_full_data_pos += full_dims[2];
        }
        data_pos += dims[1] * dims[2];
        full_data_pos += full_dims[1] * full_dims[2];
    }
}

template <class T>
void test(string filename, const vector<size_t>& dims, int retrieve_mode, double tolerance){
    size_t num_elements = 0;
    auto data_ori = MGARD::readfile<T>(filename.c_str(), num_elements);
    vector<size_t> recompose_dims(dims);
    size_t recompose_times = 0;
    auto data_recomp = test_reposition<T>(retrieve_mode, tolerance, recompose_dims, recompose_times);
    // 0 is the finest level in recomposition
    if(recompose_dims[0] != dims[0]){
        T * data_recomp_full = (T *) malloc(num_elements * sizeof(T));
        memset(data_recomp_full, 0, num_elements * sizeof(T));
        switch(dims.size()){
            case 1:
                {
                    put_back_data_1d(data_recomp, recompose_dims, data_recomp_full, dims);
                    break;
                }
            case 2:
                {
                    put_back_data_2d(data_recomp, recompose_dims, data_recomp_full, dims);
                    break;
                }
            case 3:
                {
                    put_back_data_3d(data_recomp, recompose_dims, data_recomp_full, dims);
                    break;
                }
            default:
                cerr << "Only 3 dimensions are supported in this test\n";
                exit(0);                
        }
        // recompose data to the original dimension
        MGARD::Recomposer<T> recomposer;
        int interpolation_times = (int)log2(dims[0]/recompose_dims[0]) + 1;
        recomposer.recompose_with_hierarchical_basis(data_recomp_full, dims, interpolation_times);
        MGARD::print_statistics(data_ori.data(), data_recomp_full, num_elements);
        free(data_recomp_full);
    }
    else{
        MGARD::print_statistics(data_ori.data(), data_recomp, num_elements);
    }
    free(data_recomp);
}

int main(int argc, char ** argv){
    string filename = string(argv[1]);
    int type = atoi(argv[2]); // 0 for float, 1 for double
    int retrieve_mode = atoi(argv[3]);   // 1 for max_e, 2 for squared error, 3 for PSNR
    if((retrieve_mode > 3) || (retrieve_mode < 1)) retrieve_mode = 1;
    double tolerance = atof(argv[4]);   // error tolerance 
    const int num_dims = atoi(argv[5]);
    vector<size_t> dims(num_dims);
    for(int i=0; i<dims.size(); i++){
       dims[i] = atoi(argv[6 + i]);
       cout << dims[i] << " ";
    }
    cout << endl;
    switch(type){
        case 0:
            {
                test<float>(filename, dims, retrieve_mode, tolerance);
                break;
            }
        case 1:
            {
                test<double>(filename, dims, retrieve_mode, tolerance);
                break;
            }
        default:
            cerr << "Only 0 (float) and 1 (double) are implemented in this test\n";
            exit(0);
    }
    return 0;
}