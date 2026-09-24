// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include <boost/date_time/posix_time/posix_time.hpp>

namespace s1 {
struct BurstMeta {

    boost::posix_time::ptime az_time;
    std::vector<int> first_valid_sample;
};

struct CalibMeta {
    boost::posix_time::ptime az_time;
    std::vector<float> sigma;
    std::vector<float> gamma;
    std::vector<float> beta;
    std::vector<int> pixel;
};

struct S1Metadata {
    int lines_per_burst;
    std::vector<BurstMeta> bursts;
    std::vector<CalibMeta> calib;
};
} // namespace s1