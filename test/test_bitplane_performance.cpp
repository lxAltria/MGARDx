#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include "refactor.hpp"
#include "bitstream.h"
#include <bitset>

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
    // num_elements = 8;
    // for(int i=0; i<num_elements; i++){
    //     data[i] = 1.2345;
    // }
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
    long int sum = 0;
    unsigned char * byte_encoder = buffer;
    for(int i=0; i<num_elements; i++){
        T cur_data = ldexp(data[i], prec - 1 - level_exp);
        long int fix_point = (long int) cur_data;
        bool sign = data[i] < 0;
        unsigned int fp = sign ? -fix_point : +fix_point;
        // {
        //     bitset<32> x(fp);
        //     cout << x << endl;
        // }
        int bits = 31;
        for(int j=0; j<4; j++){
            unsigned char tmp = 0;
            for(int k=7; k>=0; k--){
                tmp += ((fp >> bits) & 1) << (7-k);
                bits --;

            }
            sum += tmp;
            // *(byte_encoder ++) = tmp;
        }
        
    }
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Byteplane encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    cout << sum << endl;
    
    sum = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    unsigned char * test_pos = buffer;
    for(int i=0; i<num_elements/8; i++){
        unsigned char tmp[32] = {0};
        {
            T cur_data;
            long int fix_point;
            bool sign;

            cur_data = ldexp(data[8*i + 0], prec - 1 - level_exp);
            fix_point = (long int) cur_data;
            sign = data[8*i + 0] < 0;
            unsigned int fp0 = sign ? -fix_point : +fix_point;

            cur_data = ldexp(data[8*i + 1], prec - 1 - level_exp);
            fix_point = (long int) cur_data;
            sign = data[8*i + 1] < 0;
            unsigned int fp1 = sign ? -fix_point : +fix_point;

            cur_data = ldexp(data[8*i + 2], prec - 1 - level_exp);
            fix_point = (long int) cur_data;
            sign = data[8*i + 2] < 0;
            unsigned int fp2 = sign ? -fix_point : +fix_point;

            cur_data = ldexp(data[8*i + 3], prec - 1 - level_exp);
            fix_point = (long int) cur_data;
            sign = data[8*i + 3] < 0;
            unsigned int fp3 = sign ? -fix_point : +fix_point;

            cur_data = ldexp(data[8*i + 4], prec - 1 - level_exp);
            fix_point = (long int) cur_data;
            sign = data[8*i + 4] < 0;
            unsigned int fp4 = sign ? -fix_point : +fix_point;

            cur_data = ldexp(data[8*i + 5], prec - 1 - level_exp);
            fix_point = (long int) cur_data;
            sign = data[8*i + 5] < 0;
            unsigned int fp5 = sign ? -fix_point : +fix_point;

            cur_data = ldexp(data[8*i + 6], prec - 1 - level_exp);
            fix_point = (long int) cur_data;
            sign = data[8*i + 6] < 0;
            unsigned int fp6 = sign ? -fix_point : +fix_point;

            cur_data = ldexp(data[8*i + 7], prec - 1 - level_exp);
            fix_point = (long int) cur_data;
            sign = data[8*i + 7] < 0;
            unsigned int fp7 = sign ? -fix_point : +fix_point;

            int bits;

            bits = 31;
            for(int k=0; k<4; k++){
                for (int l=0; l<8; l++) {
                    tmp[31 - bits] += ((fp0 >> bits) & 1) << 0;
                    bits --;
                }
            }
            bits = 31;
            for(int k=0; k<4; k++){
                for (int l=0; l<8; l++) {
                    tmp[31 - bits] += ((fp1 >> bits) & 1) << 1;
                    bits --;
                }
            }
            bits = 31;
            for(int k=0; k<4; k++){
                for (int l=0; l<8; l++) {
                    tmp[31 - bits] += ((fp2 >> bits) & 1) << 2;
                    bits --;
                }
            }
            bits = 31;
            for(int k=0; k<4; k++){
                for (int l=0; l<8; l++) {
                    tmp[31 - bits] += ((fp3 >> bits) & 1) << 3;
                    bits --;
                }
            }
            bits = 31;
            for(int k=0; k<4; k++){
                for (int l=0; l<8; l++) {
                    tmp[31 - bits] += ((fp4 >> bits) & 1) << 4;
                    bits --;
                }
            }
            bits = 31;
            for(int k=0; k<4; k++){
                for (int l=0; l<8; l++) {
                    tmp[31 - bits] += ((fp5 >> bits) & 1) << 5;
                    bits --;
                }
            }
            bits = 31;
            for(int k=0; k<4; k++){
                for (int l=0; l<8; l++) {
                    tmp[31 - bits] += ((fp6 >> bits) & 1) << 6;
                    bits --;
                }
            }
            bits = 31;
            for(int k=0; k<4; k++){
                for (int l=0; l<8; l++) {
                    tmp[31 - bits] += ((fp7 >> bits) & 1) << 7;
                    bits --;
                }
            }
        }
        for(int k=0; k<32; k++){
            *(test_pos ++) = tmp[k];
        }
    }
    test_pos = buffer;
    // for(int i=0; i<4; i++){
    //     for(int j=0; j<8; j++){
    //         bitset<8> x(*test_pos++);
    //         cout << x << endl;
    //     }
    //     cout << endl;
    // }
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Byteplane encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;

    int num_level_component = 32;
    int level_component_size = num_elements / 8 + 1;
    vector<unsigned char *> byte_encoders;
    for(int i=0; i<num_level_component; i++){
        unsigned char * buffer = (unsigned char *) malloc(level_component_size);
        byte_encoders.push_back(buffer);
    }    
    err = clock_gettime(CLOCK_REALTIME, &start);
    size_t buffer_index = REFACTOR::byte_wise_direct_encoding_unrolled(buffer, num_elements, level_exp, num_level_component, byte_encoders);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Bitplane unrolling encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    for(int k=1; k<11; k++){
        // int err = clock_gettime(CLOCK_REALTIME, &start);        
        REFACTOR::RunlengthEncoder rle;
        for(int i=0; i<buffer_index; i++){
            unsigned char datum = byte_encoders[k][i];
            for(int j=0; j<8; j++){
                rle.encode(datum & 1);
                datum >>= 1;            
            }
        }
        rle.flush();
    }
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "RLE encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;

    err = clock_gettime(CLOCK_REALTIME, &start);
    int num_rle = 10;
    vector<vector<int>> counts = vector<vector<int>>(num_rle, vector<int>());
    vector<bool> last_bit = vector<bool>(num_rle, 0);
    vector<int> count = vector<int>(num_rle, 0);
    for(int i=0; i<num_elements; i++){
        T cur_data = ldexp(data[i], prec - 1 - level_exp);
        long int fix_point = (long int) cur_data;
        bool sign = data[i] < 0;
        unsigned int fp = sign ? -fix_point : +fix_point;
        // for(int bits = 30; bits >= 0; bits --){
        //     sum += ((fp >> bits) & 1);
        // }
        bool cur_bit;
        int bits = 30;
        for(int j=0; j<num_rle; j++){
            cur_bit = ((fp >> bits) & 1);
            bits --;
            if(cur_bit == last_bit[j]){
                count[j] ++;
            }
            else{
                last_bit[j] = cur_bit;
                counts[j].push_back(count[j]);
                count[j] = 0;
            }
        }
    }
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "num_elements = " << num_elements << ", counts = " << counts.size() << endl;
    // cout << *max_element(counts.begin(), counts.end()) << endl;
    cout << "RLE " << num_rle << " bitplane encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;


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