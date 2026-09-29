// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include <string>

#include "../util/proj_utils.hpp"

struct AbsTime {
    int64_t timestamp_us;
};

/*
struct RelativeTime {
    double timepoint_s;
};
*/

inline AbsTime parse_abs_time(const std::string& iso_str)
{
    // 2025-12-10T16:05:31.643322
    if (iso_str.size() != 26) {
        fmt::print("{} - parse fail\n", iso_str);
        SARTCPP_ASSERT(false);
    }
    std::tm tm { };
    long us;

    if (std::sscanf(iso_str.c_str(), "%4d-%2d-%2dT%2d:%2d:%2d.%6ld",
            &tm.tm_year, &tm.tm_mon, &tm.tm_mday,
            &tm.tm_hour, &tm.tm_min, &tm.tm_sec, &us)
        != 7) {
        fmt::print("{} - parse fail\n", iso_str);
        SARTCPP_ASSERT(false);
    }

    tm.tm_year -= 1900;
    tm.tm_mon -= 1;
    time_t ts = timegm(&tm);

    return { ts * 1000'000 + us };
}

inline std::string abstime_to_str(AbsTime abstime)
{
    time_t s = abstime.timestamp_us / 1000000;
    int64_t us = abstime.timestamp_us % 1000000;
    char buf[50];
    std::tm tm;
    gmtime_r(&s, &tm);
    strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S UTC", &tm);
    // todo iso string format... without boost
    ctime_r(&s, buf);
    std::string ret = buf;
    return ret + " + " + std::to_string(us) + " us";
}