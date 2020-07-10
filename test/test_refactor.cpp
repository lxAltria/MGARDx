#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <sys/stat.h>
#include "decompose.hpp"
#include "recompose.hpp"
#include "refactor.hpp"
#include "error_est.hpp"

using namespace std;

inline bool file_exist(const std::string& filename) {
  struct stat buffer;   
  return (stat (filename.c_str(), &buffer) == 0); 
}

template <class T>
void test_refactor(vector<T>& data, const vector<size_t>& dims, int target_level, int option, int mode, int reorganization){
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    MGARD::Decomposer<T> decomposer;
    decomposer.decompose(data.data(), dims, target_level);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Decomposition time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    // create metadata
    REFACTOR::Metadata<T> metadata(target_level);
    // set encoding option
    metadata.option = option;
    // set data reorganization
    metadata.data_reorganization = reorganization;
    // set error mode
    metadata.set_mode(mode);
    auto components = REFACTOR::level_centric_data_refactor(data.data(), target_level, dims, metadata);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Refactor time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    size_t total_size = 0;
    const vector<size_t>& level_elements = metadata.level_elements;
    const vector<T>& level_error_bounds = metadata.level_error_bounds;    
    const auto& level_errors = metadata.get_level_errors();
    unsigned char * refactored_data = NULL;
    if(metadata.data_reorganization){
        if(metadata.mode == MAX_ERROR) refactored_data = REFACTOR::refactored_data_reorganization_max_error(components, metadata.component_sizes, level_errors, metadata.order, total_size);
        else if(metadata.mode == SQUARED_ERROR) refactored_data = REFACTOR::refactored_data_reorganization_squared_error(dims.size(), components, metadata.component_sizes, level_errors, level_elements, metadata.order, total_size);
        else{
            cerr << "No such mode exist! Exit." << endl;
            exit(0);
        }
    }
    else{
        refactored_data = REFACTOR::refactored_data_reorganization_direct(components, metadata.component_sizes, metadata.order, total_size);
    }
    MGARD::writefile<unsigned char>("refactor_data/refactored.dat", refactored_data, total_size);
    // write metadata
    metadata.to_file(string("refactor_data/metadata").c_str());
    for(int i=0; i<components.size(); i++){
        for(int j=0; j<components[i].size(); j++){
            free(components[i][j]);
        }
    }
    free(refactored_data);
    // cout << "level elements: ";
    // for(int i=0; i<=target_level; i++){
    //     cout << level_elements[i] << " ";
    // }
    // cout << endl;
    // cout << "level errors: ";
    // for(int i=0; i<=target_level; i++){
    //     cout << level_error_bounds[i] << " ";
    // }
    // cout << endl << endl;
}

template <class T>
T * test_reposition(const vector<size_t>& dims, vector<size_t>& recompose_dims, size_t& recompose_times, double tolerance){
    REFACTOR::Metadata<T> metadata;
    metadata.from_file(string("refactor_data/metadata").c_str());
    int target_level = metadata.level_elements.size() - 1;
    size_t num_bytes = 0;
    const vector<vector<double>>& level_errors = metadata.get_level_errors();
    const vector<vector<size_t>>& level_sizes = metadata.component_sizes;
    const vector<int>& order = metadata.order;
    cout << "init level components\n";
    vector<int> num_intra_level_components(target_level + 1, 0);
    if(metadata.mode == SQUARED_ERROR){
        // change tolerance from psnr to se
        double psnr = tolerance;
        double value_range = 39.5582 + 53.0226;
        size_t num = 100*500*500;
        tolerance = (value_range / pow(10, psnr/20))*(value_range / pow(10, psnr/20)) * num;
        cout << "tolerance = " << tolerance << ", mean = " << tolerance / num << endl;
    }
    size_t retrieved_size = REFACTOR::interpret_reading_size(dims.size(), level_sizes, level_errors, order, metadata.mode, tolerance, num_intra_level_components);
    for(int i=0; i<=target_level; i++){
        if(num_intra_level_components[i] != 0){
            // increment number of level components because 0 and 1 are grouped as 1
            num_intra_level_components[i] ++;
        }
    }
    cout << "retrieved_size = " << retrieved_size << endl;
    unsigned char * refactored_data = MGARD::readfile_pointer<unsigned char>("refactor_data/refactored.dat", num_bytes);
    auto components = REFACTOR::read_reorganized_data(refactored_data, level_sizes, order, retrieved_size);
    int recomposed_level = 0;
    for(int i=0; i<=target_level; i++){
        if(num_intra_level_components[target_level - i] != 0){
            recomposed_level = i;
            break;
        }
    }
    cout << "Recompose to level " << recomposed_level << endl;;
    recompose_times = target_level - recomposed_level;
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    T * data = REFACTOR::level_centric_data_reposition<T>(components, metadata, target_level, recompose_times, num_intra_level_components, recompose_dims);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Reposition time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    free(refactored_data);
    err = clock_gettime(CLOCK_REALTIME, &start);
    MGARD::Recomposer<T> recomposer;
    recomposer.recompose(data, recompose_dims, recompose_times);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Recomposition time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    cout << "Recomposed dims: ";
    size_t num_elements = 1;
    for(int i=0; i<recompose_dims.size(); i++){
        cout << recompose_dims[i] << " ";
        num_elements *= recompose_dims[i];
    }
    cout << endl;
    MGARD::writefile(string("mgard.recomposed").c_str(), data, num_elements);
    num_elements = 1;
    for(const auto& d:dims){
        num_elements *= d;
    }
    cout << "Compression ratio = " << num_elements * sizeof(T) * 1.0 / retrieved_size << endl;
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
void test(string filename, const vector<size_t>& dims, int target_level, int option, int mode, double tolerance, int reorganization){
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    auto data_ori(data);
    test_refactor(data, dims, target_level, option, mode, reorganization);
    vector<size_t> recompose_dims(dims);
    size_t recompose_times = 0;
    auto data_recomp = test_reposition<T>(dims, recompose_dims, recompose_times, tolerance);
    // 0 is the finest level in recomposition
    if(recompose_times != target_level){
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
        recomposer.recompose_with_hierarchical_basis(data_recomp_full, dims, target_level - recompose_times);
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
    int target_level = atoi(argv[3]);
    int option = atoi(argv[4]); // 0 for direct, 1 for rle, 2 for hybrid
    if((option > 2) || (option < 0)) option = 0;
    int mode = atoi(argv[5]);   // 0 for max_e, 1 for squared error
    if((mode > 1) || (mode < 0)) option = 0;
    double tolerance = atof(argv[6]);   // error tolerance 
    int reorganization = atoi(argv[7]);   // enable data reorganization 
    const int num_dims = atoi(argv[8]);
    vector<size_t> dims(num_dims);
    for(int i=0; i<dims.size(); i++){
       dims[i] = atoi(argv[9 + i]);
       cout << dims[i] << " ";
    }
    cout << endl;
    switch(type){
        case 0:
            {
                test<float>(filename, dims, target_level, option, mode, tolerance, reorganization);
                break;
            }
        case 1:
            {
                test<double>(filename, dims, target_level, option, mode, tolerance, reorganization);
                break;
            }
        default:
            cerr << "Only 0 (float) and 1 (double) are implemented in this test\n";
            exit(0);
    }
    return 0;
}