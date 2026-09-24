// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include "orbit.hpp"

struct RTCParams {
    // int algo_type;
    double x_oversamp_ratio;
    double y_oversamp_ratio;
};

struct RTCDimensions {
    int dem_x_size;
    int dem_y_size;

    int az_max_idx;
    int rg_max_idx;
};

struct RTCCalcArgs {
    // double rg_spacing;

    double rg_m_to_pix;
    double az_t_to_pix;
    double slant_range_first_sample;
    double lon_start;
    double lat_start;
    double pixel_spacing_x;
    double pixel_spacing_y;
    double x_oversamp_ratio;
    double y_oversamp_ratio;
    float no_data_value;
    OrbitPolyInterpolator opi;
};
