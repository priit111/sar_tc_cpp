// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include <boost/date_time/posix_time/ptime.hpp>
#include <vector>

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
    boost::posix_time::ptime first_line_time;
    int range_size;
    int azimuth_size;
    double wavelength;
    double frequency;
    std::vector<OSV> osv;
    std::string raster_path;

    double calc_az_tp(int index) const
    {
        return line_time_interval * index;
    }

    double calc_tp_from_dt(boost::posix_time::ptime dt) const
    {
        return (first_line_time - dt).total_microseconds() * 1e-6;
    }
};
