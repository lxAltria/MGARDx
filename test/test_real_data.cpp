#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "decompose.hpp"
#include "recompose.hpp"

using namespace std;

int main(int argc, char ** argv){
	size_t num_elements = 0;
	auto data = MGARD::readfile<float>(argv[1], num_elements);
	auto data_ori(data);
    double tolerance = atof(argv[2]);
	const int target_level = atoi(argv[3]);
	// vector<size_t> dims(1, num_elements);
    const int num_dims = atoi(argv[4]);
    cout << "data size: " << num_elements << ", dims: ";
	vector<size_t> dims(num_dims);
    for(int i=0; i<dims.size(); i++){
	   dims[i] = atoi(argv[5 + i]);
       cout << dims[i] << " ";
	}
    cout << endl;
	double max_val = data[0];
	double min_val = data[0];
	double max_abs = fabs(data[0]);
	for(int i=0; i<data.size(); i++){
		if(data_ori[i] > max_val) max_val = data_ori[i];
		if(data_ori[i] < min_val) min_val = data_ori[i];
		if(fabs(data_ori[i]) > max_abs) max_abs = fabs(data_ori[i]);
	}
	double eb = tolerance;
	cout << "Required eb = " << tolerance << endl;

	MGARD::Decomposer<float> decomposer;
	size_t compressed_size = 0;
    auto compressed_data = decomposer.compress(data.data(), dims, target_level, eb, compressed_size);
    cerr << "compressed_size = " << compressed_size << endl;
    MGARD::writefile((string(argv[1]) + ".mgard").c_str(), data.data(), num_elements);
	// decomposer.decompose(data.data(), dims, target_level);
	// direct coefficient removal
	// for(int i=0; i<dims[0]; i++){
	// 	for(int j=0; j<dims[1]; j++){
	// 		if((i >= 226) || (j >= 451)) data[i * dims[1] + j] = 0;
	// 	}
	// }
	MGARD::Recomposer<float> recomposer;
	auto data_dec = recomposer.decompress(compressed_data, compressed_size, dims, target_level);
    free(compressed_data);
    // recomposer.recompose(data.data(), dims, target_level);
    // auto data_dec = data.data();
	MGARD::writefile((string(argv[1]) + ".mgard.out").c_str(), data.data(), num_elements);

	double max_err = 0;
	int pos = 0;
	double mse = 0;
	for(int i=0; i<data.size(); i++){
		double err = data_ori[i] - data_dec[i];
		mse += err * err;
		if(fabs(err) > max_err){
			pos = i;
			max_err = fabs(err);
		}
	}
	mse /= data.size();
	double psnr = 20 * log10((max_val - min_val) / sqrt(mse));
	cerr << "Max value = " << max_val << ", min value = " << min_val << endl;
	cerr << "Max error = " << max_err << ", pos = " << pos << endl;
	cerr << "MSE = " << mse << ", PSNR = " << psnr << endl;
    if(data_dec) free(data_dec);

}