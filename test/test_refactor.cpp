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

using namespace std;

inline bool file_exist(const std::string& filename) {
  struct stat buffer;   
  return (stat (filename.c_str(), &buffer) == 0); 
}

template <class T>
void test_refactor(vector<T>& data, const vector<size_t>& dims, int target_level){
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
    // whether to enable lossless compression on leading zeros
    bool with_compression = true;    
    auto components = REFACTOR::level_centric_data_refactor(data.data(), target_level, dims, metadata, with_compression);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Refactor time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    // write metadat
    // MGARD::writefile(string("refactor_data/metadata").c_str(), metadata, (target_level + 1) * (sizeof(size_t) + sizeof(T)));
    metadata.to_file(string("refactor_data/metadata").c_str());
    const vector<size_t>& level_elements = metadata.level_elements;
    const vector<T>& level_error_bounds = metadata.level_error_bounds;    
    cout << "level elements: ";
    for(int i=0; i<=target_level; i++){
        cout << level_elements[i] << " ";
    }
    cout << endl;
    cout << "level errors: ";
    for(int i=0; i<=target_level; i++){
        cout << level_error_bounds[i] << " ";
    }
    cout << endl << endl;
    for(int i=0; i<=target_level; i++){
        if(with_compression && (components[i].size() > 1)){
            size_t * bitplane_sizes = reinterpret_cast<size_t*>(components[i][0]);
            for(int j=0; j<components[i].size() - 1; j++){
                MGARD::writefile<unsigned char>(("refactor_data/level_" + to_string(target_level - i) + "_" + to_string(j)).c_str(), components[i][j+1], bitplane_sizes[j]);
                free(components[i][j+1]);
            }
            MGARD::writefile<unsigned char>(("refactor_data/level_" + to_string(target_level - i) + "_metadata").c_str(), components[i][0] + components[i].size() * sizeof(size_t), bitplane_sizes[components[i].size() - 1]);
            free(components[i][0]);
        }
        else{
            for(int j=0; j<components[i].size(); j++){
                size_t level_component_size = (level_elements[i] * sizeof(T) - 1) / components[i].size() + 1;
                MGARD::writefile<unsigned char>(("refactor_data/level_" + to_string(target_level - i) + "_" + to_string(j)).c_str(), components[i][j], level_component_size);
                free(components[i][j]);
            }
        }
    }
}

template <class T>
T * test_reposition(const vector<size_t>& dims, int target_recompose_level, vector<size_t>& recompose_dims){
    vector<vector<unsigned char*>> components;
    // auto metadata = MGARD::readfile_pointer<unsigned char>(string("refactor_data/metadata").c_str(), tmp_size);
    REFACTOR::Metadata<T> metadata;
    metadata.from_file(string("refactor_data/metadata").c_str());
    int target_level = metadata.level_elements.size() - 1;
    // auto intra_level_components = MGARD::readfile_pointer<unsigned char>(string("refactor_data/metadata").c_str(), tmp_size);
    vector<int> num_intra_level_components(target_level + 1, 32);
    num_intra_level_components[0] = 1;
    // num_intra_level_components[target_level] = 16;
    target_recompose_level = target_level - target_recompose_level;
    vector<bool> with_compression(target_recompose_level + 1, false);
    for(int i=0; i<=target_recompose_level; i++){
        size_t tmp_size = 0;
        // check whether "metadata" exists
        with_compression[i] = file_exist(("refactor_data/level_" + to_string(target_level - i) + "_metadata").c_str());
        vector<unsigned char*> level_component;
        if(with_compression[i]){
            // push back metadata
            auto intra_level_metadata = MGARD::readfile_pointer<unsigned char>(("refactor_data/level_" + to_string(target_level - i) + "_metadata").c_str(), tmp_size);
            level_component.push_back(intra_level_metadata);            
        }
        for(int j=0; j<num_intra_level_components[i]; j++){
            // push back components
            auto intra_level_component = MGARD::readfile_pointer<unsigned char>(("refactor_data/level_" + to_string(target_level - i) + "_" + to_string(j)).c_str(), tmp_size);
            level_component.push_back(intra_level_component);
        }
        components.push_back(level_component);
    }
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    T * data = REFACTOR::level_centric_data_reposition<T>(components, metadata, target_level, target_recompose_level, num_intra_level_components, recompose_dims, with_compression);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Reposition time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    for(int i=0; i<components.size(); i++){
        for(int j=0; j<components[i].size(); j++){
            free(components[i][j]);
        }
    }
    err = clock_gettime(CLOCK_REALTIME, &start);
    MGARD::Recomposer<T> recomposer;
    recomposer.recompose(data, recompose_dims, target_recompose_level);
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
void test(string filename, const vector<size_t>& dims, int target_level, int target_recompose_level){
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    auto data_ori(data);
    test_refactor(data, dims, target_level);
    vector<size_t> recompose_dims(dims);
    auto data_recomp = test_reposition<T>(dims, target_recompose_level, recompose_dims);
    // 0 is the finest level in recomposition
    if(target_recompose_level != 0){
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
        recomposer.recompose_with_hierarchical_basis(data_recomp_full, dims, target_recompose_level);
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
    int target_recompose_level = atoi(argv[4]);
    if(target_level < target_recompose_level) target_recompose_level = target_level;
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
                test<float>(filename, dims, target_level, target_recompose_level);
                break;
            }
        case 1:
            {
                test<double>(filename, dims, target_level, target_recompose_level);
                break;
            }
        default:
            cerr << "Only 0 (float) and 1 (double) are implemented in this test\n";
            exit(0);
    }
    return 0;
}