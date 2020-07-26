#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <sys/stat.h>
#include "refactor.hpp"
#include "error_est.hpp"

using namespace std;

inline bool file_exist(const std::string& filename) {
  struct stat buffer;   
  return (stat (filename.c_str(), &buffer) == 0); 
}

template <class T>
void test_refactor(string filename, const vector<size_t>& dims, int target_level, int option, int reorganization){
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    unsigned char * refactored_data = NULL;
    REFACTOR::Metadata<T> metadata = REFACTOR::multigrid_data_refactor(data, dims, target_level, option, reorganization, &refactored_data);
    MGARD::writefile<unsigned char>("refactor_data/refactored.dat", refactored_data, metadata.total_encoded_size);
    metadata.to_file(string("refactor_data/metadata").c_str());
    cout << "data written to refactor_data/refactored.dat" << endl;
    cout << "metadata written to refactor_data/metadata" << endl;
    free(refactored_data);
}

int main(int argc, char ** argv){
    string filename = string(argv[1]);
    int type = atoi(argv[2]); // 0 for float, 1 for double
    int target_level = atoi(argv[3]);
    int option = atoi(argv[4]); // 0 for direct, 1 for direct+sign postpone, 2 for rle, 3 for hybrid, 4 for hybrid+sign postpone 
    if((option > 4) || (option < 0)) option = 0;
    int reorganization = atoi(argv[5]);   // enable data reorganization 
    const int num_dims = atoi(argv[6]);
    vector<size_t> dims(num_dims);
    for(int i=0; i<dims.size(); i++){
       dims[i] = atoi(argv[7 + i]);
       cout << dims[i] << " ";
    }
    cout << endl;
    switch(type){
        case 0:
            {
                test_refactor<float>(filename, dims, target_level, option, reorganization);
                break;
            }
        case 1:
            {
                test_refactor<double>(filename, dims, target_level, option, reorganization);
                break;
            }
        default:
            cerr << "Only 0 (float) and 1 (double) are implemented in this test\n";
            exit(0);
    }
    return 0;
}