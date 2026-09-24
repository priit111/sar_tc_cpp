// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include <fmt/format.h>
//#include <fmt/ranges.h>

#include "../sar/sar_metadata.hpp"
#include "../util/memory_raster.hpp"
#include "s1_metadata.hpp"

namespace s1 {

inline std::vector<size_t> deburst_find_az_offsets(const SARMetadata& sar_meta, const S1Metadata& s1_meta)
{
    const double lti = sar_meta.line_time_interval;

    auto first_time = s1_meta.bursts.front().az_time;
    int burst_sz = s1_meta.lines_per_burst;

    std::vector<double> burst_starts;
    std::vector<double> burst_ends;
    std::vector<double> burst_first_line;
    for (int burst_idx = 0; burst_idx < s1_meta.bursts.size(); burst_idx++) {

        const BurstMeta& bm = s1_meta.bursts[burst_idx];
        size_t valid_begin = burst_sz;
        size_t valid_end = burst_sz;
        for (size_t i = 0; i < bm.first_valid_sample.size(); i++) {
            if (valid_begin == burst_sz && bm.first_valid_sample[i] != -1) {
                valid_begin = i;
            }
            if (valid_begin != burst_sz && bm.first_valid_sample[i] == -1) {
                valid_end = i - 1;
                break;
            }
        }

        double burst_start_time = (bm.az_time - first_time).total_microseconds() * 1e-6;
        // fmt::print("burst = {} {}\n", valid_begin, valid_end);
        double start_time = burst_start_time + valid_begin * lti;
        double end_time = burst_start_time + valid_end * lti;
        if (burst_idx > 0) {
            double mid_tp = (burst_ends[burst_idx - 1] + burst_start_time) / 2.0;
            start_time = mid_tp;
            burst_ends[burst_idx - 1] = mid_tp;
        }
        burst_first_line.push_back(burst_start_time);
        burst_starts.push_back(start_time);
        burst_ends.push_back(end_time);
    }

    std::vector<size_t> offsets;
    {
        int burst_idx = 0;
        // fmt::print("BEGs = {}\n", burst_starts);
        // fmt::print("ENDS = {}\n", burst_ends);
        burst_starts[0] = 0;
        int cnt = 0;
        while (burst_idx != burst_starts.size()) {
            double tp = cnt * lti;
            if (tp > burst_ends[burst_idx]) {
                burst_idx++;
                continue;
            }
            double dt = tp - burst_first_line[burst_idx];
            int idx = std::round(dt / lti);
            idx += burst_sz * burst_idx;
            offsets.push_back(idx);
            cnt++;
        }
    }
    return offsets;
}

inline void deburst(const SARMetadata& sar_meta, const S1Metadata& s1_meta, MemoryRaster<float>& in_out)
{

    /*
    const double lti = sar_meta.line_time_interval;

    auto first_time = s1_meta.bursts.front().az_time;
    int burst_sz = s1_meta.lines_per_burst;

    std::vector<double> burst_starts;
    std::vector<double> burst_ends;
    std::vector<double> burst_first_line;
    for (int burst_idx = 0; burst_idx < s1_meta.bursts.size(); burst_idx++) {

        const BurstMeta& bm = s1_meta.bursts[burst_idx];
        size_t valid_begin = burst_sz;
        size_t valid_end = burst_sz;
        for (size_t i = 0; i < bm.first_valid_sample.size(); i++) {
            if (valid_begin == burst_sz && bm.first_valid_sample[i] != -1) {
                valid_begin = i;
            }
            if (valid_begin != burst_sz && bm.first_valid_sample[i] == -1) {
                valid_end = i - 1;
                break;
            }
        }

        double burst_start_time = (bm.az_time - first_time).total_microseconds() * 1e-6;
        //fmt::print("burst = {} {}\n", valid_begin, valid_end);
        double start_time = burst_start_time + valid_begin * lti;
        double end_time = burst_start_time + valid_end * lti;
        if (burst_idx > 0) {
            double mid_tp = (burst_ends[burst_idx - 1] + burst_start_time) / 2.0;
            start_time = mid_tp;
            burst_ends[burst_idx - 1] = mid_tp;
        }
        burst_first_line.push_back(burst_start_time);
        burst_starts.push_back(start_time);
        burst_ends.push_back(end_time);
    }

    std::vector<size_t> offsets;
    {
        int burst_idx = 0;
        //fmt::print("BEGs = {}\n", burst_starts);
        //fmt::print("ENDS = {}\n", burst_ends);
        burst_starts[0] = 0;
        int cnt = 0;
        while (burst_idx != burst_starts.size()) {
            double tp = cnt * lti;
            if (tp > burst_ends[burst_idx]) {
                burst_idx++;
                continue;
            }
            double dt = tp - burst_first_line[burst_idx];
            int idx = std::round(dt /  lti);
            idx += burst_sz * burst_idx;
            offsets.push_back(idx);
            cnt++;
        }
    }*/

    const std::vector<size_t> offsets = deburst_find_az_offsets(sar_meta, s1_meta);

    const size_t x_size = in_out.m_x_size;
    const size_t y_size = offsets.size();

    float* in_p = in_out.m_data;

    for (size_t y = 0; y < y_size; y++) {
        size_t y_offset = offsets[y];
        if (y == y_offset) {
            continue;
        }
        size_t in_idx = y_offset * x_size;
        size_t out_idx = y * x_size;
        memcpy(&in_p[out_idx], &in_p[in_idx], x_size * sizeof(float));
    }

    in_out.resize_y(offsets.size());
}
} // namespace s1