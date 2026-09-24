// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#include "multilook.hpp"

#include "../proj_conf.hpp"

void multilook(SARMetadata& sar_meta, MemoryRaster<float>& memory_raster, int range_multilook)
{

    TimeBlock tb("range multilook\n");
    int y_size = memory_raster.m_y_size;
    int x_size = memory_raster.m_x_size;

    int new_x_size = x_size / range_multilook;

    fmt::print("size = {} -> {}\n", x_size, new_x_size);
    int out_idx = 0;
    for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < new_x_size; x++) {

            float sum = 0.0f;
            const float* data_in = &memory_raster.m_data[y * x_size + x * range_multilook];
            for (int i = 0; i < range_multilook; i++) {
                sum += data_in[i];
            }

            sum /= range_multilook;

            memory_raster.m_data[out_idx] = sum;
            out_idx++;
        }
    }

    memory_raster.resize_xy(new_x_size, y_size);

    sar_meta.range_spacing *= range_multilook;
    sar_meta.range_size = new_x_size;
}
