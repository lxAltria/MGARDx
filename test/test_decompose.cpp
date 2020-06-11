#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include "decompose.hpp"
#include "recompose.hpp"

using namespace std;

template <class T>
void test_decompose(vector<T>& data, const vector<size_t>& dims, int target_level){
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    MGARD::Decomposer<T> decomposer;
    decomposer.decompose(data.data(), dims, target_level);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Decomposition time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
}

template <class T>
void test_recompose(vector<T>& data, const vector<size_t>& dims, int target_level){
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    MGARD::Recomposer<T> recomposer;
    recomposer.recompose(data.data(), dims, target_level);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Recomposition time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
}

template <class T>
void test(string filename, const vector<size_t>& dims, int target_level){
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    auto data_ori(data);
    test_decompose(data, dims, target_level);
    test_recompose(data, dims, target_level);
    MGARD::print_statistics(data_ori.data(), data.data(), num_elements);
}

int main(int argc, char ** argv){
    string filename = string(argv[1]);
    int type = atoi(argv[2]); // 0 for float, 1 for double
    int target_level = atoi(argv[3]);
    const int num_dims = atoi(argv[4]);
    vector<size_t> dims(num_dims);
    for(int i=0; i<dims.size(); i++){
       dims[i] = atoi(argv[5 + i]);
       cout << dims[i] << " ";
    }
    cout << endl;
    switch(type){
        case 0:
            {
                test<float>(filename, dims, target_level);
                break;
            }
        case 1:
            {
                test<double>(filename, dims, target_level);
                break;
            }
        default:
            cerr << "Only 0 (float) and 1 (double) are implemented in this test\n";
            exit(0);
    }
    return 0;
}