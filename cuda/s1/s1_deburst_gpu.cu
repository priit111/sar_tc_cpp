//
// Created by priit on 9/24/26.
//

#include "s1_deburst_gpu.hpp"

#include "../../s1/s1_deburst.hpp"

namespace s1 {

void deburst_gpu(const SARMetadata& sar_meta, const S1Metadata& s1_meta, DeviceRaster<float>& in_out)
{
    TimeBlock t("s1 deburst gpu");
    std::vector<size_t> offsets = deburst_find_az_offsets(sar_meta, s1_meta);

    const size_t x_size = in_out.m_x_size;
    const size_t y_size = offsets.size();

    float* in_p = in_out.m_d_data;

    for (size_t y = 0; y < y_size; y++) {
        size_t y_offset = offsets[y];
        if (y == y_offset) {
            continue;
        }
        size_t in_idx = y_offset * x_size;
        size_t out_idx = y * x_size;
        cudaMemcpy(&in_p[out_idx], &in_p[in_idx], x_size * sizeof(float), cudaMemcpyDeviceToDevice);
    }

    in_out.resize_y(offsets.size());
}

}