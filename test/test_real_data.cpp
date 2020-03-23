#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <fstream>
#include "decompose.hpp"
#include "recompose.hpp"

using namespace std;

template<typename Type>
std::vector<Type> readfile(const char *file, size_t &num) {
    std::ifstream fin(file, std::ios::binary);
    if (!fin) {
        std::cout << " Error, Couldn't find the file" << "\n";
        return std::vector<Type>();
    }
    fin.seekg(0, std::ios::end);
    const size_t num_elements = fin.tellg() / sizeof(Type);
    fin.seekg(0, std::ios::beg);
    auto data = std::vector<Type>(num_elements);
    fin.read(reinterpret_cast<char *>(&data[0]), num_elements * sizeof(Type));
    fin.close();
    num = num_elements;
    return data;
}

template<typename Type>
void writefile(const char *file, Type *data, size_t num_elements) {
    std::ofstream fout(file, std::ios::binary);
    fout.write(reinterpret_cast<const char *>(&data[0]), num_elements * sizeof(Type));
    fout.close();
}

int main(int argc, char ** argv){
	size_t num_elements = 0;
	auto data = readfile<float>(argv[1] , num_elements);
	const int target_level = atoi(argv[2]);
	cerr << "Origin data: " << endl;
	for(int i=0; i<20; i++){
		cerr << data[i] << " ";
	}
	cerr << endl;
	MGARD::Decomposer<float> decomposer;
	vector<size_t> dims(1, num_elements);
	decomposer.decompose(data.data(), dims, target_level);

	MGARD::Recomposer<float> recomposer;
	recomposer.recompose(data.data(), dims, target_level);
	cerr << "Recomposed data: " << endl;
	for(int i=0; i<20; i++){
		cerr << data[i] << " ";
	}
	cerr << endl;
}