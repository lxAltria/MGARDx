#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include "utils.hpp"
#include "reorder.hpp"
#include "misc.hpp"
#include "correction.hpp"

using namespace MGARD;

template<class T>
std::vector<double> compute_operator_norm_3D(T * Q, const vector<size_t>& dims, size_t target_level, double s){
    const int default_batch_size = 32;
    size_t n1 = dims[0];
    size_t n2 = dims[1];
    size_t n3 = dims[2];
    size_t dim0_stride = dims[1] * dims[2];
    size_t dim1_stride = dims[2];
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n = n1 * n2 * n3;
    const T h = 1;
    T * data_buffer = (T *) malloc(n * sizeof(T));
    T * Q_copy = (T *) malloc(n * sizeof(T));
    T * load_v_buffer = (T *) malloc(default_batch_size * (*max_element(dims.begin(), dims.end())) * sizeof(T));
    vector<double> operater_norm(target_level + 1, 0);
    // compute eta_L
    data_reorder_3D(Q, data_buffer, n1, n2, n3, dim0_stride, dim1_stride);
    data_copy_3D(Q_copy, Q, n1, n2, n3, dim0_stride, dim1_stride);
    compute_correction_3D(Q, data_buffer, load_v_buffer, n1, n2, n3, n1_nodal, h, dim0_stride, dim1_stride, default_batch_size);
    double prev_level_norm = dot_product_3D(Q, Q_copy, n1, n2, n3, dim0_stride, dim1_stride);
    for(int i=0; i<target_level; i++){
        restriction_3D(Q, n1, n2, n3, dim0_stride, dim1_stride);
        data_reorder_3D(Q, data_buffer, n1, n2, n3, dim0_stride, dim1_stride);
        data_copy_3D(Q_copy, Q, n1, n2, n3, dim0_stride, dim1_stride);    
        compute_correction_3D(Q, data_buffer, load_v_buffer, n1, n2, n3, n1_nodal, h, dim0_stride, dim1_stride, default_batch_size);
        n1 = (n1 >> 1) + 1;
        n2 = (n2 >> 1) + 1;
        n3 = (n3 >> 1) + 1;
        double level_norm = dot_product_3D(Q, Q_copy, n1, n2, n3, dim0_stride, dim1_stride);
        operater_norm[target_level - i] = pow(2, -2*s*(target_level - i)) * (prev_level_norm - level_norm);
        prev_level_norm = level_norm;
    }
    operater_norm[0] = prev_level_norm;
    free(data_buffer);
    free(Q_copy);
    free(load_v_buffer);
    return operater_norm;
}


template <class T>
void test_operator_norm(vector<T>& data, const vector<size_t>& dims, int target_level, double s){
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    auto norm = compute_operator_norm_3D(data.data(), dims, target_level, s);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Computing operator norm time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    cout << "Level operator norms:" << endl;
    for(const auto& norm_l:norm){
        cout << norm_l << endl;
    }
    cout << endl;
}

template <class T>
void test(string filename, const vector<size_t>& dims, int target_level, double s){
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    test_operator_norm(data, dims, target_level, s);
}

int main(int argc, char ** argv){
    string filename = string(argv[1]);
    int type = atoi(argv[2]); // 0 for float, 1 for double
    int target_level = atoi(argv[3]);
    double s = atof(argv[4]);
    const int num_dims = atoi(argv[5]);
    vector<size_t> dims(num_dims);
    for(int i=0; i<dims.size(); i++){
       dims[i] = atoi(argv[6 + i]);
       cout << dims[i] << " ";
    }
    cout << endl;
    fflush(stdout);
    switch(type){
        case 0:
            {
                test<float>(filename, dims, target_level, s);
                break;
            }
        case 1:
            {
                test<double>(filename, dims, target_level, s);
                break;
            }
        default:
            cerr << "Only 0 (float) and 1 (double) are implemented in this test\n";
            exit(0);
    }
    return 0;
}