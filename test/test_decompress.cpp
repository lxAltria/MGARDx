#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include "recompose.hpp"

using namespace std;

template <class T>
T * test_decompress(const unsigned char * compressed, size_t compressed_size, const vector<size_t>& dims){
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    MGARD::Recomposer<T> recomposer;
    auto data_dec = recomposer.decompress(compressed, compressed_size, dims);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cerr << "Decompression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    return data_dec;
}

template <class T>
void test(string ori_filename, string compressed_filename, const vector<size_t>& dims){
    size_t num_elements = 0;
    size_t compressed_size = 0;
    auto compressed = MGARD::readfile<unsigned char>(compressed_filename.c_str(), compressed_size);
    auto data_ori = MGARD::readfile<float>(ori_filename.c_str(), num_elements);
    auto data_dec = test_decompress<float>(compressed.data(), compressed_size, dims);
    MGARD::print_statistics(data_ori.data(), data_dec, num_elements, compressed_size);
    MGARD::writefile((compressed_filename + ".out").c_str(), data_dec, num_elements);
    free(data_dec);    
}

int main(int argc, char ** argv){
    string ori_filename = string(argv[1]);
    string compressed_filename = string(argv[2]);
    int type = atoi(argv[3]); // 0 for float, 1 for double
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
                test<float>(ori_filename, compressed_filename, dims);
                break;
            }
        case 1:
            {
                test<double>(ori_filename, compressed_filename, dims);
                break;
            }
        default:
            cerr << "Only 0 (float) and 1 (double) are implemented in this test\n";
            exit(0);
    }
    return 0;
}