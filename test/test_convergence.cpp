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
void test_decompose(vector<T>& data, const vector<size_t>& dims, int target_level, int block_size){
	// 3D data by default
	assert(dims.size() == 3);
    vector<size_t> strides(3, 0);
    {
		size_t stride = 1;
		for(int i=dims.size()-1; i>=0; i--){
			strides[i] = stride;
			stride *= dims[i];
		}    	
    }
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    int nx = (dims[0] - 1) / block_size + 1;
    int ny = (dims[1] - 1) / block_size + 1;
    int nz = (dims[2] - 1) / block_size + 1;
    T * x_data_pos = data.data();
    vector<size_t> block_dims(3, 0);
    for(int i=0; i<nx; i++){
    	T * y_data_pos = x_data_pos;
    	block_dims[0] = (i < nx - 1) ? block_size : (dims[0] - i*block_size);
    	for(int j=0; j<ny; j++){
    		T * z_data_pos = y_data_pos;
	    	block_dims[1] = (j < ny - 1) ? block_size : (dims[1] - j*block_size);
    		for(int k=0; k<nz; k++){
		    	block_dims[2] = (k < nz - 1) ? block_size : (dims[2] - k*block_size);
    			cout << "dimensions = " << block_dims[0] << " " << block_dims[1] << " " << block_dims[2] << endl;
			    MGARD::Decomposer<T> decomposer;
		    	decomposer.decompose(z_data_pos, block_dims, target_level, false, strides);
		    	z_data_pos += block_size;
    		}
    		y_data_pos += block_size * dims[2];
    	}
    	x_data_pos += block_size * dims[1] * dims[2];
    }
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Decomposition time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
}

template <class T>
void test_recompose(vector<T>& data, const vector<size_t>& dims, int target_level, int block_size){
    vector<size_t> strides(3, 0);
    {
		size_t stride = 1;
		for(int i=dims.size()-1; i>=0; i--){
			strides[i] = stride;
			stride *= dims[i];
		}    	
    }
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    int nx = (dims[0] - 1) / block_size + 1;
    int ny = (dims[1] - 1) / block_size + 1;
    int nz = (dims[2] - 1) / block_size + 1;
    T * x_data_pos = data.data();
    vector<size_t> block_dims(3, 0);
    for(int i=0; i<nx; i++){
    	T * y_data_pos = x_data_pos;
    	block_dims[0] = (i < nx - 1) ? block_size : (dims[0] - i*block_size);
    	for(int j=0; j<ny; j++){
    		T * z_data_pos = y_data_pos;
	    	block_dims[1] = (j < ny - 1) ? block_size : (dims[1] - j*block_size);
    		for(int k=0; k<nz; k++){
		    	block_dims[2] = (k < nz - 1) ? block_size : (dims[2] - k*block_size);
			    MGARD::Recomposer<T> recomposer;
		    	recomposer.recompose(z_data_pos, block_dims, target_level, false, strides);
		    	z_data_pos += block_size;
    		}
    		y_data_pos += block_size * dims[2];
    	}
    	x_data_pos += block_size * dims[1] * dims[2];
    }
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Recomposition time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
}

template <class T>
void verify(const T * data_ori, const T * data_dec, const vector<size_t>& dims, int block_size, T isovalve){
    int nx = (dims[0] - 1) / block_size + 1;
    int ny = (dims[1] - 1) / block_size + 1;
    int nz = (dims[2] - 1) / block_size + 1;
    const T * x_data_pos = data_ori;
    vector<size_t> block_dims(3, 0);
    for(int i=0; i<nx; i++){
    	const T * y_data_pos = x_data_pos;
    	block_dims[0] = (i < nx - 1) ? block_size : (dims[0] - i*block_size);
    	for(int j=0; j<ny; j++){
    		const T * z_data_pos = y_data_pos;
	    	block_dims[1] = (j < ny - 1) ? block_size : (dims[1] - j*block_size);
    		for(int k=0; k<nz; k++){
		    	block_dims[2] = (k < nz - 1) ? block_size : (dims[2] - k*block_size);
		    	// verify
		    	T max_err = 0;
		    	double l2_error = 0;
			    int iso_change_count = 0;
		    	const T * xx_data_pos = z_data_pos;
		    	for(int ii=0; ii<block_dims[0]; ii++){
		    		const T * yy_data_pos = xx_data_pos;
		    		for(int jj=0; jj<block_dims[1]; jj++){
		    			const T * zz_data_pos = yy_data_pos;
		    			for(int kk=0; kk<block_dims[2]; kk++){
		    				T val = *zz_data_pos;
		    				ptrdiff_t offset = zz_data_pos - data_ori;
		    				T dec_val = *(data_dec + offset);
		    				T diff = fabs(val - dec_val);
		    				if(diff > max_err) max_err = diff;
		    				l2_error += diff * diff;
		    				if((val - isovalve)*(dec_val - isovalve) < 0){
		    					iso_change_count ++;
		    				}
		    				zz_data_pos ++;
		    			}
		    			yy_data_pos += dims[2];
		    		}
		    		xx_data_pos += dims[1] * dims[2];
		    	}
				cout << "i, j, k = " << i << " " << j << " " << k << endl;
				cout << "dimensions = " << block_dims[0] << " " << block_dims[1] << " " << block_dims[2] << endl;
				cout << "max error = " << max_err << endl;
				cout << "rmse = " << sqrt(l2_error / (block_dims[0]*block_dims[1]*block_dims[2])) << endl;
				cout << "iso_change_count = " << iso_change_count << endl;
		    	z_data_pos += block_size;
    		}
    		y_data_pos += block_size * dims[2];
    	}
    	x_data_pos += block_size * dims[1] * dims[2];
    }
}

int recompose_levels[200] = {
		4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 0, 0, 0, 1, 0, 4, 4,
       4, 4, 0, 0, 0, 0, 4, 0, 4, 4, 4, 4, 0, 0, 0, 0, 1, 0, 4, 4, 4, 4,
       0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 0, 0, 0, 1, 1, 0, 4, 4, 4, 4, 4, 0,
       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
       4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
       4, 4, 4, 4, 4, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
       4, 4, 4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 0, 4, 4, 4,
       4, 4, 1, 4, 0, 1, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 0, 4, 4, 4, 4, 4,
       4, 4
};
template <class T>
void clear_data(vector<T>& data, const vector<size_t>& dims, int target_level, int total_level, int block_size){
    int nx = (dims[0] - 1) / block_size + 1;
    int ny = (dims[1] - 1) / block_size + 1;
    int nz = (dims[2] - 1) / block_size + 1;
    T * x_data_pos = data.data();
    vector<size_t> block_dims(3, 0);
    for(int i=0; i<nx; i++){
    	T * y_data_pos = x_data_pos;
    	block_dims[0] = (i < nx - 1) ? block_size : (dims[0] - i*block_size);
    	for(int j=0; j<ny; j++){
    		T * z_data_pos = y_data_pos;
	    	block_dims[1] = (j < ny - 1) ? block_size : (dims[1] - j*block_size);
    		for(int k=0; k<nz; k++){
		    	block_dims[2] = (k < nz - 1) ? block_size : (dims[2] - k*block_size);
			    // clear data in each block: reset fine-level data to 0
			    // specific decomposition
			    // int target_level = recompose_levels[i*ny*nz + j*nz + k] + 1;
			    // int target_level = total_level - 1;
			    auto level_dims = MGARD::init_levels(block_dims, total_level);
			    auto current_dims = level_dims[target_level];
			    cout << " current_dims = " << current_dims[0] << " " << current_dims[1] << " " << current_dims[2] << endl;
			    T * xx_data_pos = z_data_pos;
			    for(int ii=0; ii<block_dims[0]; ii++){
			    	T * yy_data_pos = xx_data_pos;
			    	for(int jj=0; jj<block_dims[1]; jj++){
			    		T * zz_data_pos = yy_data_pos;
			    		for(int kk=0; kk<block_dims[2]; kk++){
			    			if((ii < current_dims[0]) && (jj < current_dims[1]) && (kk < current_dims[2])){
			    				// keep low level data
			    			}
			    			else{
			    				*zz_data_pos = 0;
			    			}
			    			zz_data_pos ++;
			    		}
			    		yy_data_pos += dims[2];
			    	}
			    	xx_data_pos += dims[1] * dims[2];
			    }
		    	z_data_pos += block_size;
    		}
    		y_data_pos += block_size * dims[2];
    	}
    	x_data_pos += block_size * dims[1] * dims[2];
    }
}

template <class T>
void test(string filename, const vector<size_t>& dims, int target_level, int block_size, int option=0){
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    auto data_ori(data);
	T isovalve = 0;
	{
		T max_val = data_ori[0];
		T min_val = data_ori[0];
		for(int i=1; i<dims[0]*dims[1]*dims[2]; i++){
			if(data_ori[i] > max_val) max_val = data_ori[i];
			if(data_ori[i] < min_val) min_val = data_ori[i]; 
		}
		isovalve = (min_val + max_val) / 2;
		printf("isovalue = %.14f\n", isovalve);
	}
    test_decompose(data, dims, target_level, block_size);
 //    {
	//     // check specific recomposition
	//     clear_data(data, dims, target_level, target_level, block_size);
	//     test_recompose(data, dims, target_level, block_size);
	//     verify(data_ori.data(), data.data(), dims, block_size, isovalve);
	//     MGARD::writefile("recomposed_multiresolution.dat", data.data(), data.size());
	// }
    vector<T> data_pre;
    for(int i=0; i<=target_level; i++){
    	cout << "Recompose level = " << i << endl;
	    auto data_cleaned(data);
	    clear_data(data_cleaned, dims, i, target_level, block_size);
	    test_recompose(data_cleaned, dims, target_level, block_size);
	    string filename = "recomposed_level_" + to_string(i); 
	    MGARD::writefile(filename.c_str(), data_cleaned.data(), data_cleaned.size());
	    if(option == 0){
	    	// evaluate actual error
		    verify(data_ori.data(), data_cleaned.data(), dims, block_size, isovalve);    	
	    }
	    else{
	    	// evaluate convergence error
		    if(i > 0){
		    	verify(data_cleaned.data(), data_pre.data(), dims, block_size, isovalve);
		    }
		    data_pre = data_cleaned;	    	
	    }
    }
}

int main(int argc, char ** argv){
	int arg_pos = 1;
    string filename = string(argv[arg_pos++]);
    // floating-point data by default
    int target_level = atoi(argv[arg_pos++]);
    const int num_dims = atoi(argv[arg_pos++]);
    vector<size_t> dims(num_dims);
    for(int i=0; i<dims.size(); i++){
       dims[i] = atoi(argv[arg_pos++]);
       cout << dims[i] << " ";
    }
    cout << endl;
    int block_size = atoi(argv[arg_pos++]);
    int option = atoi(argv[arg_pos++]);
    cout << "option = " << option << endl;
    test<float>(filename, dims, target_level, block_size, option);
    return 0;
}