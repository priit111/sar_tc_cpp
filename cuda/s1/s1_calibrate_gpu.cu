//
// Created by priit on 9/24/26.
//

#include "s1_calibrate_gpu.hpp"

__global__ void BetaCalibrate(IQ16orF32* in_out_data, int x_size, int y_size, float beta_coeff)
{
    const int x = threadIdx.x + blockDim.x * blockIdx.x;
    const int y = threadIdx.y + blockDim.y * blockIdx.y;
    const int data_idx = y * x_size + x;

    if (x < x_size && y < y_size) {
        IQ16 iq = in_out_data[data_idx].iq16;

        float I = iq.i;
        float Q = iq.q;

        in_out_data[data_idx].f32 = (I * I + Q * Q) * beta_coeff;
    }
}

namespace s1 {
void calibrate_gpu(const S1Metadata& s1_meta, DeviceRaster<IQ16orF32>& data_in_out)
{

    TimeBlock t("s1 calibrate gpu");
    float beta = s1_meta.calib.front().beta.front(); // an hack for now, assume beta = constant

    int x_size = data_in_out.m_x_size;
    int y_size = data_in_out.m_y_size;
    dim3 block_size(16, 16);
    dim3 grid_size((x_size + 15) / 16, (y_size + 15) / 16);

    BetaCalibrate<<<grid_size, block_size>>>(data_in_out.m_d_data, x_size, y_size, 1.0f / (beta * beta));
    CUDA_VERIFY_CHECK(cudaDeviceSynchronize()); // not needed at this point, but helps with debugging
    CUDA_VERIFY_CHECK(cudaGetLastError());
}
}