#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include "decompose.hpp"
#include "recompose.hpp"
#include "utils.hpp"

using namespace std;

template <class T>
vector<T> rand_vec(int n, T scale=1){
    vector<T> result(n);
    for(int i=0; i<n; i++){
        result[i] = scale * ( (rand() * 2.0 / RAND_MAX) - 1);
    }
    return result;
}

int main(int argc, char ** argv){
    const size_t height = atoi(argv[1]);
    const size_t depth = atoi(argv[2]);
    const size_t width = atoi(argv[3]);
    int target_level = atoi(argv[4]);
    vector<size_t> dims{height, depth, width};
    const int n = height * depth * width;
    cout << n << endl;
    vector<double> data = rand_vec<double>(n);
    for(int i=0; i<n; i++){
        data[i] = i;
    }
    auto data_(data);
    // for(int i=0; i<height; i++){
    //     for(int j=0; j<depth; j++){
    //         for(int k=0; k<width; k++){
    //             cout << data[i * depth * width + j * width + k] << " ";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }
    // cout << endl << endl;
    MGARD::Decomposer<double> decomposer;
    decomposer.decompose(data.data(), dims, target_level);

    // for(int i=0; i<height; i++){
    //     for(int j=0; j<depth; j++){
    //         for(int k=0; k<width; k++){
    //             cout << data[i * depth * width + j * width + k] << " ";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }
    // cout << endl << endl;

    MGARD::Recomposer<double> recomposer;
    recomposer.recompose(data.data(), dims, target_level);

    // for(int i=0; i<height; i++){
    //     for(int j=0; j<depth; j++){
    //         for(int k=0; k<width; k++){
    //             cout << data[i * depth * width + j * width + k] << " ";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }
    // cout << endl << endl;

    double error = 0;
    for(int i=0; i<n; i++){
        auto err = fabs(data_[i] - data[i]);
        if(err > error) error = err;
    }
    cout << error << endl;

}