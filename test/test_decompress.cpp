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
void print_statistics(const T * data_ori, const T * data_dec, size_t data_size, size_t compressed_size){
    double max_val = data_ori[0];
    double min_val = data_ori[0];
    double max_abs = fabs(data_ori[0]);
    for(int i=0; i<data_size; i++){
        if(data_ori[i] > max_val) max_val = data_ori[i];
        if(data_ori[i] < min_val) min_val = data_ori[i];
        if(fabs(data_ori[i]) > max_abs) max_abs = fabs(data_ori[i]);
    }
    double max_err = 0;
    int pos = 0;
    double mse = 0;
    for(int i=0; i<data_size; i++){
        double err = data_ori[i] - data_dec[i];
        mse += err * err;
        if(fabs(err) > max_err){
            pos = i;
            max_err = fabs(err);
        }
    }
    mse /= data_size;
    double psnr = 20 * log10((max_val - min_val) / sqrt(mse));
    cout << "Max value = " << max_val << ", min value = " << min_val << endl;
    cout << "Max error = " << max_err << ", pos = " << pos << endl;
    cout << "MSE = " << mse << ", PSNR = " << psnr << endl;
    cout << "Compression ratio = " << data_size * sizeof(T) * 1.0 / compressed_size << endl;
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
    size_t num_elements = 0;
    size_t compressed_size = 0;
    auto compressed = MGARD::readfile<unsigned char>(compressed_filename.c_str(), compressed_size);
    switch(type){
        case 0:
            {
                auto data_ori = MGARD::readfile<float>(ori_filename.c_str(), num_elements);
                auto data_dec = test_decompress<float>(compressed.data(), compressed_size, dims);
                print_statistics(data_ori.data(), data_dec, num_elements, compressed_size);
                MGARD::writefile((compressed_filename + ".out").c_str(), data_dec, num_elements);
                free(data_dec);
                break;
            }
        case 1:
            {
                auto data_ori = MGARD::readfile<double>(ori_filename.c_str(), num_elements);
                auto data_dec = test_decompress<double>(compressed.data(), compressed_size, dims);
                print_statistics(data_ori.data(), data_dec, num_elements, compressed_size);
                MGARD::writefile((compressed_filename + ".out").c_str(), data_dec, num_elements);
                free(data_dec);
                break;
            }
        default:
            cerr << "Only 0 (float) and 1 (double) are implemented in this test\n";
            exit(0);
    }
    return 0;
}