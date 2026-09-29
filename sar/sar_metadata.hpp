// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include <vector>

#include "../util/proj_utils.hpp"
#include "time_type.hpp"

struct OSV {
    double tp;
    double xp;
    double yp;
    double zp;
    double xv;
    double yv;
    double zv;
};

inline OSV InterpolateOrbit(const OSV* osv, int n, double tp)
{

    OSV r = { };
    r.tp = tp;
    for (int i = 0; i < n; i++) {
        double mult = 1;
        for (int j = 0; j < n; j++) {
            if (i == j)
                continue;

            double xj = osv[j].tp;
            double xi = osv[i].tp;
            mult *= (tp - xj) / (xi - xj);
        }

        r.xp += mult * osv[i].xp;
        r.yp += mult * osv[i].yp;
        r.zp += mult * osv[i].zp;
        r.xv += mult * osv[i].xv;
        r.yv += mult * osv[i].yv;
        r.zv += mult * osv[i].zv;
    }

    return r;
}

struct SARMetadata {
    double azimuth_spacing;
    double range_spacing;
    double slant_range_first_sample;
    double line_time_interval;
    AbsTime first_line_time;
    int range_size;
    int azimuth_size;
    double wavelength;
    double frequency;
    std::vector<OSV> osv;
    std::string raster_path;

    double incidence_angle_begin;
    double incidence_angle_end;

    double calc_az_tp(int index) const
    {
        return line_time_interval * index;
    }

    void cut_y(int new_size, AbsTime new_first_line_time)
    {
        double orig_last_tp = line_time_interval * azimuth_size;
        double delta = (new_first_line_time.timestamp_us - first_line_time.timestamp_us) * 1e-6;
        SARTCPP_ASSERT(delta >= 0);
        double new_last_tp = delta + new_size * line_time_interval;
        SARTCPP_ASSERT(new_last_tp <= orig_last_tp);

        first_line_time = new_first_line_time;

        for (auto& osv_e : osv) {
            osv_e.tp -= delta;
        }
    }
};
