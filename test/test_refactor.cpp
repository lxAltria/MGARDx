#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include "decompose.hpp"
#include "recompose.hpp"
#include "refactor.hpp"

using namespace std;

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
    auto components = MGARD::level_centric_data_refactor(data.data(), target_level, dims);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Refactor time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    MGARD::writefile(string("refactor.metadata").c_str(), components[0], (target_level + 1) * (sizeof(size_t) + sizeof(T)));
    size_t * level_elements = reinterpret_cast<size_t*>(components[0]);
    T * level_error_bounds = reinterpret_cast<T*>(components[0] + (target_level + 1) * sizeof(size_t));    
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
        MGARD::writefile<unsigned char>(("refactor.level_" + to_string(target_level - i)).c_str(), components[i + 1], level_elements[i] * sizeof(T));
        free(components[i + 1]);
    }
    free(components[0]);
}
template <class T>
T * test_reposition(const vector<size_t>& dims, int target_recompose_level){
    vector<unsigned char*> components;
    size_t tmp_size = 0;
    auto metadata = MGARD::readfile_pointer<unsigned char>(string("refactor.metadata").c_str(), tmp_size);
    int target_level = tmp_size / (sizeof(T) + sizeof(size_t)) - 1;
    cout << target_level << endl;
    size_t * level_elements = reinterpret_cast<size_t*>(metadata);
    components.push_back(metadata);
    target_recompose_level = target_level - target_recompose_level;
    for(int i=0; i<=target_recompose_level; i++){
        auto level_component = MGARD::readfile_pointer<unsigned char>(("refactor.level_" + to_string(target_level - i)).c_str(), tmp_size);
        components.push_back(level_component);
    }
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    auto recompose_dims(dims);
    T * data = MGARD::level_centric_data_reposition<T>(components, target_level, target_recompose_level, recompose_dims);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Reposition time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    for(int i=0; i<components.size(); i++){
        free(components[i]);
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
    size_t num_elements = 0;
    switch(type){
        case 0:
            {
                auto data = MGARD::readfile<float>(filename.c_str(), num_elements);
                test_refactor(data, dims, target_level);
                auto data_recomp = test_reposition<float>(dims, target_recompose_level);
                free(data_recomp);
                break;
            }
        case 1:
            {
                auto data = MGARD::readfile<double>(filename.c_str(), num_elements);
                test_refactor(data, dims, target_level);
                auto data_recomp = test_reposition<double>(dims, target_recompose_level);
                free(data_recomp);
                break;
            }
        default:
            cerr << "Only 0 (float) and 1 (double) are implemented in this test\n";
            exit(0);
    }
    return 0;
}