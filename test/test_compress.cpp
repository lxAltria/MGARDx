#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include "utils.hpp"
#include "decompose.hpp"
#include "recompose.hpp"

using namespace std;

template <class T>
unsigned char * test_compress(vector<T>& data, const vector<size_t>& dims, int target_level, double eb, size_t& compressed_size){
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    MGARD::Decomposer<T> decomposer;
    auto compressed_data = decomposer.compress(data.data(), dims, target_level, eb, compressed_size);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Compression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    return compressed_data;
}

int main(int argc, char ** argv){
    string filename = string(argv[1]);
    int type = atoi(argv[2]); // 0 for float, 1 for double
    double tolerance = atof(argv[3]);
    const int target_level = atoi(argv[4]);
    const int num_dims = atoi(argv[5]);
    vector<size_t> dims(num_dims);
    for(int i=0; i<dims.size(); i++){
       dims[i] = atoi(argv[6 + i]);
       cout << dims[i] << " ";
    }
    cout << endl;
    double eb = tolerance;
    cout << "Required eb = " << tolerance << endl;
    size_t num_elements = 0;
    size_t compressed_size = 0;
    unsigned char * compressed = NULL;
    switch(type){
        case 0:
            {
                auto data = MGARD::readfile<float>(filename.c_str(), num_elements);
                compressed = test_compress(data, dims, target_level, eb, compressed_size);
                break;
            }
        case 1:
            {
                auto data = MGARD::readfile<double>(filename.c_str(), num_elements);
                compressed = test_compress(data, dims, target_level, eb, compressed_size);
                break;
            }
        default:
            cerr << "Only 0 (float) and 1 (double) are implemented in this test\n";
            exit(0);
    }
    MGARD::writefile((filename + ".mgard").c_str(), compressed, compressed_size);
    free(compressed);
    return 0;
}