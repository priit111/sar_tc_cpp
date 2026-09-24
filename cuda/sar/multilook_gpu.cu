// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#include "multilook_gpu.hpp"


__global__ void multilook_kernel(const float* data_in, float* data_out,int in_x_size, int out_x_size, int out_y_size, int range_multilook)
{
    const int x = threadIdx.x + blockDim.x * blockIdx.x;
    const int y = threadIdx.y + blockDim.y * blockIdx.y;

    if (x > out_x_size || y > out_y_size) {
        return;
    }

    float sum = 0.0f;
    const float* ptr = &data_in[y * in_x_size + range_multilook * x];
    for (int i = 0; i < range_multilook; i++) {
        sum += ptr[i];
    }
    sum /= range_multilook;

    data_out[y * out_x_size + x] = sum;
}



void multilook_gpu(SARMetadata& sar_meta, DeviceRaster<float>& memory_raster, int range_multilook)
{
    TimeBlock tb("range gpu multilook\n");
    int in_y_size = memory_raster.m_y_size;
    int in_x_size = memory_raster.m_x_size;


    int out_x_size = in_x_size / range_multilook;
    int out_y_size = in_y_size;

    fmt::print("size = {} -> {}\n", in_x_size, out_x_size);

    DeviceRaster<float> multilooked_buf;
    multilooked_buf.init(out_x_size, out_y_size);

    dim3 block_size(16, 16);
    dim3 grid_size((out_x_size + 15) / 16, (out_y_size + 15) / 16);

    multilook_kernel<<<grid_size, block_size>>>(memory_raster.m_d_data, multilooked_buf.m_d_data, in_x_size, out_x_size, out_y_size, range_multilook);

    sar_meta.range_spacing *= range_multilook;
    sar_meta.range_size = out_x_size;;

    memory_raster = std::move(multilooked_buf);

    CUDA_VERIFY_CHECK(cudaDeviceSynchronize()); // not needed at this point, but helps with debugging
    CUDA_VERIFY_CHECK(cudaGetLastError());
}