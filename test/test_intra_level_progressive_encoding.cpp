#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include "refactor.hpp"

using namespace std;

template <class T>
void test(string filename, int recompose_level_intra){
    struct timespec start, end;
    int err = 0;
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    T max_val = 0;
    for(const auto& d:data){
        if(fabs(d) > max_val) max_val = fabs(d);
    }
    int level_exp = 0;
    frexp(max_val, &level_exp);
    cout << max_val << " " << level_exp << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    auto level_components = MGARD::progressive_encoding(data.data(), num_elements, level_exp);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Progressive encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    T * data_recomp = MGARD::progressive_decoding<T>(level_components, num_elements, level_exp, recompose_level_intra);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Progressive decoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    MGARD::print_statistics(data.data(), data_recomp, num_elements);
    free(data_recomp);
    for(int i=0; i<level_components.size(); i++){
        free(level_components[i]);
    }
}

int main(int argc, char **argv){
    string filename = string(argv[1]);
    int recompose_level_intra = atoi(argv[2]);
    test<float>(filename, recompose_level_intra);
}