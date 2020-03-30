#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "decompose.hpp"
#include "recompose.hpp"
#include <quantizer/Quantizer.hpp>
#include <encoder/HuffmanEncoder.hpp>
#include "zstd.h"

using namespace std;

unsigned long sz_lossless_compress(unsigned char *data, size_t dataLength) {
    unsigned long outSize = 0;
    size_t estimatedCompressedSize = 0;
    if (dataLength < 100)
        estimatedCompressedSize = 200;
    else
        estimatedCompressedSize = dataLength * 1.2;
    unsigned char * buffer = (unsigned char *) malloc(estimatedCompressedSize);
    outSize = ZSTD_compress(buffer, estimatedCompressedSize, data, dataLength,
                            3); //default setting of level is 3
    free(buffer);
    return outSize;
}

int main(int argc, char ** argv){
	size_t num_elements = 0;
	auto data = MGARD::readfile<float>(argv[1], num_elements);
	auto data_ori(data);
	const int target_level = atoi(argv[2]);
	// vector<size_t> dims(1, num_elements);
	vector<size_t> dims(2);
	dims[0] = atoi(argv[3]);
	dims[1] = atoi(argv[4]);

	double max_val = data[0];
	double min_val = data[0];
	double max_abs = fabs(data[0]);
	for(int i=0; i<data.size(); i++){
		if(data_ori[i] > max_val) max_val = data_ori[i];
		if(data_ori[i] < min_val) min_val = data_ori[i];
		if(fabs(data_ori[i]) > max_abs) max_abs = fabs(data_ori[i]);
	}
	double eb = atof(argv[5]);
	double C2 = 1 + pow(sqrt(3)/2, dims.size());
	double t = eb * (max_val - min_val) / (C2 * (target_level + 1));
	cout << "Required eb = " << eb << ", tolerance = " << t << endl;

	MGARD::Decomposer<float> decomposer;
	decomposer.decompose(data.data(), dims, target_level);

	// Huffman encoding
	{
		unsigned char * compressed = (unsigned char *) malloc(num_elements * sizeof(float));
		unsigned char * compressed_data_pos = compressed;
		auto quantizer = SZ::LinearQuantizer<float>(t);
		auto encoder = SZ::HuffmanEncoder<int>();
		vector<int> quant_inds(num_elements);
		for(int i=0; i<num_elements; i++){
			quant_inds[i] = quantizer.quantize(data[i], 0);
			data[i] = quantizer.recover(0, quant_inds[i]);
		}
		encoder.preprocess_encode(quant_inds, 4*quantizer.get_radius());
		encoder.save(compressed_data_pos);
		encoder.encode(quant_inds, compressed_data_pos);
		encoder.postprocess_encode();
		size_t compressed_length = compressed_data_pos - compressed;
		cout << "Huffman encoded length = " << compressed_length << endl;
		cout << "Compression ratio = " << num_elements * sizeof(float) * 1.0 / compressed_length << endl;
		auto final_length = sz_lossless_compress(compressed, compressed_length);
		cout << "Final length = " << final_length << endl;
		cout << "Compression ratio = " << num_elements * sizeof(float) * 1.0 / final_length << endl;
	}

	// direct coefficient removal
	// for(int i=0; i<dims[0]; i++){
	// 	for(int j=0; j<dims[1]; j++){
	// 		if((i >= 226) || (j >= 451)) data[i * dims[1] + j] = 0;
	// 	}
	// }

	MGARD::writefile((string(argv[1]) + ".mgard").c_str(), data.data(), num_elements);

	MGARD::Recomposer<float> recomposer;
	recomposer.recompose(data.data(), dims, target_level);

	MGARD::writefile((string(argv[1]) + ".mgard.out").c_str(), data.data(), num_elements);

	double max_err = 0;
	int pos = 0;
	double mse = 0;
	for(int i=0; i<data.size(); i++){
		double err = data_ori[i] - data[i];
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
	// MGARD::Recomposer<float> recomposer;
	// recomposer.recompose(data.data(), dims, target_level);
	// cerr << "Recomposed data: " << endl;
	// for(int i=0; i<20; i++){
	// 	cerr << data[i] << " ";
	// }
	// cerr << endl;
}