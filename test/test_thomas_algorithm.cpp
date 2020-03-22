#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
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
vector<T> solve_tridiagonal(vector<T>& load_v_buffer, T h, size_t n){
	// Thomas algorithm for solving M_l x = load_v
	// forward pass
	// simplified algorithm
	// b[0] = h/3, b[1:n] = h/4	
	T * d = load_v_buffer.data();
	vector<T> b(n, h/3);
	T c = h/6;
	for(int i=1; i<n; i++){
		auto w = c / b[i-1];
		b[i] = b[i] - w * c;
		d[i] = d[i] - w * d[i-1];
	}
	// backward pass
	vector<T> result(n);
	result[n-1] = d[n-1] / b[n-1];
	for(int i=n-2; i>0; i--){
		result[i] = (d[i] - c * result[i+1]) / b[i];
	}
	result[0] = (d[0] - c * result[1]) / b[0];
	return result;
}

int main(int argc, char ** argv){
	int n = atoi(argv[1]);
	double h = atof(argv[2]);
	auto d = rand_vec<double>(n);
	auto d_prime = vector<double>(d);
	// for(int i=0; i<n; i++){
	// 	cout << d[i] << " ";
	// }
	// cout << endl;
	auto x = solve_tridiagonal(d_prime, h, n);
	auto diff = (h/3 * x[0] + h/6 * x[1]) - d[0];
	cout << diff << " ";
	for(int i=1; i<n-1; i++){
		auto diff = (h/6 * x[i-1] + h/3 * x[i] + h/6 * x[i+1]) - d[i];
		cout << diff << " ";
	}
	diff = (h/6 * x[n-2] + h/3 * x[n-1]) - d[n-1];
	cout << diff << endl;
}