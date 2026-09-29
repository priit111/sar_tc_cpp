// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once


#include "../sar/time_type.hpp"
#include <map>

namespace s1 {
struct BurstMeta {

    AbsTime az_time;
    std::vector<int> first_valid_sample;
};

struct CalibMeta {
    AbsTime az_time;
    std::vector<float> sigma;
    std::vector<float> gamma;
    std::vector<float> beta;
    std::vector<int> pixel;
};


struct GeoLocationGridPoint {
    std::string az_time;
    double slrt;
    int line;
    int pixel;
    double lat;
    double lon;
    double height;
    double incidence;
    double elevation;
};

struct S1Metadata {
    int lines_per_burst;
    std::vector<BurstMeta> bursts;
    std::vector<CalibMeta> calib;

    std::map<int, std::vector<GeoLocationGridPoint>> geogrid_points;
};
} // namespace s1