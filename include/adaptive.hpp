#include <cmath>

namespace MGARD{

using namespace std;

template <class T>
bool switch_to_lorenzo(const T * data, size_t n1, size_t n2, size_t n3, size_t dim0_stride, size_t dim1_stride, double eb){
    const T * data_pos = data + dim0_stride + dim1_stride + 1;
    const int block_size = 3;
    const int skip_count = 5; // sample every 5 3x3 blocks => overhead = 1/125
    const size_t stride = block_size * skip_count;
    if((n1 < stride) || (n2 < stride) || (n3 < stride)) return true;
    size_t num_x = (n1 - 1) / stride;
    size_t num_y = (n2 - 1) / stride;
    size_t num_z = (n3 - 1) / stride;
    double lorenzo_err = 0;
    double lorenzo_noise = eb * 1.22;
    // TODO: consult Ben for detail?
    double nodal_noise_2 = eb * 0.6665/2;
    double nodal_noise_4 = eb * 0.9332/4;
    double nodal_noise_8 = eb * 1.31234/8;
    double interpolation_err = 0;
    const T * x_data_pos = data_pos;
    const vector<size_t> coefficient_stride{
        1, 2*dim1_stride + 1, 2*dim0_stride + 1, 2*dim0_stride + 2*dim1_stride + 1, // nodal_nodal_coeff
        dim1_stride, dim1_stride + 2, 2*dim0_stride + dim1_stride, 2*dim0_stride + dim1_stride + 2, // nodal_coeff_nodal
        dim0_stride, dim0_stride + 2, dim0_stride + 2*dim1_stride, dim0_stride + 2*dim1_stride + 2, // coeff_nodal_nodal
        dim1_stride + 1, 2*dim0_stride + dim1_stride + 1, // nodal_coeff_coeff
        dim0_stride + 1, dim0_stride + 2*dim1_stride + 1, // coeff_nodal_coeff
        dim0_stride + dim1_stride, dim0_stride + dim1_stride + 2, // coeff_coeff_nodal
        dim0_stride + dim1_stride + 1 // coeff_coeff_coeff
    };
    for(int i=0; i<num_x; i++){
        const T * y_data_pos = x_data_pos;
        for(int j=0; j<num_y; j++){
            const T * z_data_pos = y_data_pos;
            for(int k=0; k<num_z; k++){
                // iterate through 3^3 - 2^3 = 19 coefficients
                {
                    // nodal_nodal_coeff
                    for(int i=0; i<4; i++){
                        ptrdiff_t offset = coefficient_stride[i];
                        interpolation_err += fabs((z_data_pos[offset - 1] + z_data_pos[offset + 1]) / 2 - z_data_pos[offset]);
                        interpolation_err += nodal_noise_2;
                    }
                    // nodal_coeff_nodal
                    for(int i=4; i<8; i++){
                        ptrdiff_t offset = coefficient_stride[i];
                        interpolation_err += fabs((z_data_pos[offset - dim1_stride] + z_data_pos[offset + dim1_stride]) / 2 - z_data_pos[offset]);
                        interpolation_err += nodal_noise_2;
                    }
                    // coeff_nodal_nodal
                    for(int i=8; i<12; i++){
                        ptrdiff_t offset = coefficient_stride[i];
                        interpolation_err += fabs((z_data_pos[offset - dim0_stride] + z_data_pos[offset + dim0_stride]) / 2 - z_data_pos[offset]);
                        interpolation_err += nodal_noise_2;
                    }
                    // nodal_coeff_coeff
                    for(int i=12; i<14; i++){
                        ptrdiff_t offset = coefficient_stride[i];
                        interpolation_err += fabs((z_data_pos[offset - dim1_stride - 1] + z_data_pos[offset - dim1_stride + 1] + z_data_pos[offset + dim1_stride - 1] + z_data_pos[offset + dim1_stride + 1]) / 4 - z_data_pos[offset]);
                        interpolation_err += nodal_noise_4;
                    }
                    // coeff_nodal_coeff
                    for(int i=14; i<16; i++){
                        ptrdiff_t offset = coefficient_stride[i];
                        interpolation_err += fabs((z_data_pos[offset - dim0_stride - 1] + z_data_pos[offset - dim0_stride + 1] + z_data_pos[offset + dim0_stride - 1] + z_data_pos[offset + dim0_stride + 1]) / 4 - z_data_pos[offset]);
                        interpolation_err += nodal_noise_4;
                    }
                    // coeff_coeff_nodal
                    for(int i=16; i<18; i++){
                        ptrdiff_t offset = coefficient_stride[i];
                        interpolation_err += fabs((z_data_pos[offset - dim0_stride - dim1_stride] + z_data_pos[offset - dim0_stride + dim1_stride] + z_data_pos[offset + dim0_stride - dim1_stride] + z_data_pos[offset + dim0_stride + dim1_stride]) / 4 - z_data_pos[offset]);
                        interpolation_err += nodal_noise_4;
                    }
                    // coeff_coeff_coeff
                    for(int i=18; i<19; i++){
                        ptrdiff_t offset = coefficient_stride[i];
                        interpolation_err += fabs((z_data_pos[offset - dim0_stride - dim1_stride - 1] + z_data_pos[offset - dim0_stride - dim1_stride + 1] +
                                              z_data_pos[offset - dim0_stride + dim1_stride - 1] + z_data_pos[offset - dim0_stride + dim1_stride + 1] +
                                              z_data_pos[offset + dim0_stride - dim1_stride - 1] + z_data_pos[offset + dim0_stride - dim1_stride + 1] +
                                              z_data_pos[offset + dim0_stride + dim1_stride - 1] + z_data_pos[offset + dim0_stride + dim1_stride + 1]
                                                ) / 8 - z_data_pos[offset]);
                        interpolation_err += nodal_noise_8;
                    }
                    // lorenzo
                    for(int i=0; i<19; i++){
                        size_t offset = coefficient_stride[i];
                        T predicted = z_data_pos[offset - 1] + z_data_pos[offset - dim0_stride] + z_data_pos[offset - dim1_stride]
                                        - z_data_pos[offset - dim1_stride - 1] - z_data_pos[offset - dim0_stride - 1] - z_data_pos[offset - dim0_stride - dim1_stride]
                                        + z_data_pos[offset - dim0_stride - dim1_stride - 1];
                        lorenzo_err += fabs(predicted - z_data_pos[offset]);
                    }
                    lorenzo_err += 19 * lorenzo_noise;
                }
                z_data_pos += stride;
            }
            y_data_pos += stride * dim1_stride;
        }
        x_data_pos += stride * dim0_stride;
    }
    cout << "lorenzo_err = " << lorenzo_err << endl;
    cout << "interpolation_err = " << interpolation_err << endl;
    return lorenzo_err < interpolation_err;
}    

}