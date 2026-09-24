#pragma once

#include "../util/cuda_util.hpp"
#include <array>

struct DeviceDEM {
    DeviceBuffer<float> d_data;
    int x_size;
    int y_size;
    float no_data_value;
    std::array<double, 6> gt;

    void init(const float* h_data_in, size_t x_size, size_t y_size, float no_data_value, std::array<double, 6> gt)
    {
        this->x_size = x_size;
        this->y_size = y_size;
        this->no_data_value = no_data_value;
        this->gt = gt;

        d_data.init(x_size * y_size);
        h2d_cpy(d_data.data(), h_data_in, x_size * y_size);
    }
};
