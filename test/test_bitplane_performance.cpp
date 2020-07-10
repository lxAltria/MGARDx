#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include "refactor.hpp"
#include "bitstream.h"

using namespace std;

// A class to write data bit by bit
// Currently does not consider the case when encoding size is larger than capacity
class BitEncoder{
protected:
    bitstream* stream = NULL;  
public:
    BitEncoder(unsigned char * array, size_t size){
        stream = stream_open(array, size);
    }
    void close(){
        stream_close(stream);
    }
    void encode(bool b){
        stream_write_bit(stream, b);
    }
    void flush(){
        stream_flush(stream);
    }
    size_t size(){ return stream_size(stream); }
};

class ByteEncoder{
private:
    unsigned char * start = NULL;
    unsigned char * current = NULL;
public:
    ByteEncoder(unsigned char * s){
        start = s;
        current = s;
    }
    void encode(unsigned char c){
        *(current ++) = c;
    }
    void flush(){}
    void close(){}
    size_t size(){
        return current - start;
    }
};

template <class T>
void test(string filename){
    struct timespec start, end;
    int err = 0;
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    err = clock_gettime(CLOCK_REALTIME, &start);
    auto max_val = REFACTOR::record_level_max_value(data.data(), num_elements);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Compute max_val time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    int level_exp = 0;
    frexp(max_val, &level_exp);
    cout << max_val << " " << level_exp << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    // encode
    unsigned char * buffer = (unsigned char *) malloc(num_elements*sizeof(T));
    const int prec = 32;
    BitEncoder encoder(buffer, num_elements*sizeof(T));
    for(int i=0; i<num_elements; i++){
        T cur_data = ldexp(data[i], prec - 1 - level_exp);
        long int fix_point = (long int) cur_data;
        bool sign = data[i] < 0;
        encoder.encode(sign);
        unsigned int fp = sign ? -fix_point : +fix_point;
        for(int j=prec - 1; j>0; j--){
            encoder.encode(fp & 1);
            fp >>= 1;
        }  
    }
    encoder.flush();
    encoder.close();
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Bitplane encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;

    err = clock_gettime(CLOCK_REALTIME, &start);
    unsigned char * byte_encoder = buffer;
    for(int i=0; i<num_elements; i++){
        T cur_data = ldexp(data[8], prec - 1 - level_exp);
        long int fix_point = (long int) cur_data;
        bool sign = data[i] < 0;
        unsigned int fp = sign ? -fix_point : +fix_point;
        int bits = 31;
        for(int j=0; j<4; j++){
            unsigned char tmp = 0;
            for(int k=7; k>=0; k--){
                tmp += ((fp >> bits) & 1) << (7-k);
                bits --;
            }
            *(byte_encoder ++) = tmp;
        }
    }
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Byteplane encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    
    err = clock_gettime(CLOCK_REALTIME, &start);
    ByteEncoder b_encoder(buffer);
    for(int i=0; i<num_elements/8; i++){
        unsigned char tmp[32] = {0};
        for(int j=0; j<8; j++){
            T cur_data = ldexp(data[8*i + j], prec - 1 - level_exp);
            long int fix_point = (long int) cur_data;
            bool sign = data[i] < 0;
            unsigned int fp = sign ? -fix_point : +fix_point;
            for(int k=31; k>=0; k--){
                tmp[k] += ((fp >> k) & 1) << j;
            }
        }
        for(int k=0; k<32; k++)
            b_encoder.encode(tmp[k]);
    }
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Byteplane encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;

    free(buffer);

    err = clock_gettime(CLOCK_REALTIME, &start);
    REFACTOR::record_level_mse(data.data(), num_elements, 24, level_exp);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Compute A_l time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    // MGARD::print_statistics(data.data(), data_recomp, num_elements);
}

int main(int argc, char **argv){
    string filename = string(argv[1]);
    int recompose_level_intra = atoi(argv[2]);
    test<float>(filename);
}