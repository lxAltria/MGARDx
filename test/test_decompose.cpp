#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include "decompose.hpp"
#include "recompose.hpp"

using namespace std;

template <class T>
vector<T> rand_vec(int n, T scale=1){
	vector<T> result(n);
	for(int i=0; i<n; i++){
		result[i] = scale * ( (rand() * 2.0 / RAND_MAX) - 1);
	}
	return result;
}

template <class T>
void compute_interpolant_difference(vector<T>& data, const vector<T>& data_, int n, int target_stride){
	data.push_back(data.back());
	double interpolant = 0;
	double mse = 0;
	for(int i=0; i<n; i++){
		int residue = i % target_stride;
		int ind = i / target_stride;
		if(residue){
			double lambda = residue * 1.0 / target_stride;
			interpolant = (1 - lambda) * data[ind * target_stride] + lambda * data[(ind+1)*target_stride];
		}
		else interpolant = data[i];
		mse += (data_[i] - interpolant) * (data_[i] - interpolant);
		cerr << setprecision(4) << data_[i] - interpolant << " ";
	}
	cerr << endl;
	cerr << "MSE = " << mse << endl;
	data.pop_back();
}

int main(int argc, char ** argv){
	const int n = atoi(argv[1]);
	const int target_level = atoi(argv[2]);
	const int target_stride = 1 << target_level;
	vector<double> data(n);
	for(int i=0; i<n; i++){
		data[i] = 0.05 * (i - 5) * (i - 10) * (i - 15); 
	}
	auto data_(data);
	for(int i=0; i<n; i+=target_stride){
		cout << data[i] << " ";
	}
	cout << endl;
	// direct interpolant
	compute_interpolant_difference(data, data_, n, target_stride);
	// MGARD::decompose(data.data(), n, target_level);
	MGARD::Decomposer<double> decomposer;
	vector<size_t> dims(1, n);
	decomposer.decompose(data.data(), dims, target_level);
	// MGARD interpolant
	vector<double> data_reordered = vector<double>(n, 0);
	for(int i=0, j=0; i<n; i+=target_stride, j++){
		cout << data[j] << " ";
		data_reordered[i] = data[j];
	}
	cout << endl;
	compute_interpolant_difference(data_reordered, data_, n, target_stride);
	MGARD::Recomposer<double> recomposer;
	recomposer.recompose(data.data(), dims, target_level);
	cerr << "Origin data: " << endl;
	for(int i=0; i<n; i++){
		cerr << data_[i] << " ";
	}
	cerr << endl;
	cerr << "Recomposed data: " << endl;
	for(int i=0; i<n; i++){
		cerr << data[i] << " ";
	}
	cerr << endl;
}