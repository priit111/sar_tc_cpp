// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp


#pragma once

#include <cmath>

#include "../proj_conf.hpp"
#include "../util/math_utils.hpp"
#include "orbit.hpp"

#ifndef GEOCODE_MODE
#error "GEOCODE_MODE..."
#endif

struct TCArgs {
    // geo lat/lon
    double lon_start;
    double lat_start;
    double pixel_spacing_y;
    double pixel_spacing_x;

    float no_data_value;

    // geometry of the acquisition
    double line_time_interval;
    double range_spacing;
    double slant_range_first_sample;

    // radar
    double wavelength;

#if GEOCODE_MODE == 0
    // OSV
    const OSV* osv;
    int n_osv;

    // LUTs for the zero Doppler binary search
    const Pos3D* sat_pos_az;
    const Vec3D* sat_vel_az;
#elif GEOCODE_MODE == 1
    OrbitPolyInterpolator opi;
#endif
};

struct TCDimensions {
    int in_x_size;
    int in_y_size;
    int out_x_size;
    int out_y_size;
};

namespace {
namespace WGS84 { // NOLINT
    constexpr double A = 6378137.0; // m
    constexpr double B = 6356752.3142451794975639665996337; // 6356752.31424518; // m
    constexpr double FLAT_EARTH_COEF = 1.0 / ((A - B) / A); // 298.257223563;
    constexpr double E2 = 2.0 / FLAT_EARTH_COEF - 1.0 / (FLAT_EARTH_COEF * FLAT_EARTH_COEF);
    constexpr double EP2 = E2 / (1 - E2);
} // namespace WGS84

constexpr double DTOR = M_PI / 180.0;
constexpr double RTOD = 180 / M_PI;
}

HD_FUNC
inline Pos3D Geo2xyzWgs84(double latitude, double longitude, double altitude)
{
    double const lat = latitude * DTOR;
    double const lon = longitude * DTOR;

    double sin_lat;
    double cos_lat;
    sincos(lat, &sin_lat, &cos_lat);

    double const sinLat = sin_lat;

    double const N = (WGS84::A / sqrt(1.0 - WGS84::E2 * sinLat * sinLat));
    double const NcosLat = (N + altitude) * cos_lat;

    double sin_lon;
    double cos_lon;
    sincos(lon, &sin_lon, &cos_lon);

    Pos3D r = { };
    r.x = NcosLat * cos_lon;
    r.y = NcosLat * sin_lon;
    r.z = (N + altitude - WGS84::E2 * N) * sinLat;
    return r;
}

constexpr double INVALID_RDR_IDX = -1234567e7; // not nan due to -Ofast

#if GEOCODE_MODE == 0

HD_FUNC
inline double CalcDopplerFrequency(Pos3D earth_point, Pos3D sensor_pos, Vec3D sensor_vel, double wavelength)
{
    const Vec3D dpos = earth_point - sensor_pos;
    const auto dis = length(dpos);

    return 2.0 * dot(dpos, sensor_vel) / (dis * wavelength);
}

HD_FUNC
inline double GetEarthPointZeroDopplerTime(double line_time_interval, double wavelength, Pos3D earth_point, int n_azimuth,
    const Pos3D* sensor_position, const Vec3D* sensor_velocity)
{
    // binary search is used in finding the zero doppler time
    int lower_bound = 0;
    int upper_bound = n_azimuth - 1;
    auto lower_bound_freq = CalcDopplerFrequency(earth_point, sensor_position[lower_bound], sensor_velocity[lower_bound], wavelength);
    auto upper_bound_freq = CalcDopplerFrequency(earth_point, sensor_position[upper_bound], sensor_velocity[upper_bound], wavelength);

    if (std::abs(lower_bound_freq) < 1.0) {
        return lower_bound * line_time_interval;
    } else if (std::abs(upper_bound_freq) < 1.0) {
        return upper_bound * line_time_interval;
    } else if (lower_bound_freq * upper_bound_freq > 0.0) {
        return INVALID_RDR_IDX;
    }

    /*if constexpr (ZERO_DOPPLER_INTERCEPT) {
        double k = (upper_bound_freq - lower_bound_freq) / ((n_azimuth)*line_time_interval);
        double C = upper_bound_freq - k * n_azimuth * line_time_interval;
        double az_idx = -C / k;
        int az_idx_i = az_idx / line_time_interval;
        upper_bound = az_idx_i + 5;
        lower_bound = az_idx_i - 5;
    }*/

    // start binary search
    double mid_freq;
    while (upper_bound - lower_bound > 1) {
        const auto mid = (int)((lower_bound + upper_bound) / 2.0);
        mid_freq = sensor_velocity[mid].x * (earth_point.x - sensor_position[mid].x) + sensor_velocity[mid].y * (earth_point.y - sensor_position[mid].y) + sensor_velocity[mid].z * (earth_point.z - sensor_position[mid].z);

        if (mid_freq * lower_bound_freq > 0.0) {
            lower_bound = mid;
            lower_bound_freq = mid_freq;
        } else if (mid_freq * upper_bound_freq > 0.0) {
            upper_bound = mid;
            upper_bound_freq = mid_freq;
        } else if (mid_freq == 0.0) {
            return mid * line_time_interval;
        }
    }

    const auto y0 = lower_bound - lower_bound_freq * (upper_bound - lower_bound) / (upper_bound_freq - lower_bound_freq);

    // fmt::print(" bin search = {} , slope = {}, delta {}\n", y0 * line_time_interval, az_calc_idx, az_calc_idx - y0 * line_time_interval);
    return y0 * line_time_interval;
}

HD_FUNC
inline Pos3D GetPosition(double time, const OSV* vectors, int n_osv)
{
    int i0 = 0;
    int iN = n_osv - 1;
    if (iN > 8) {
        double t_first = vectors[0].tp;
        i0 = std::max<int>(0, (time - t_first) - 4);
        iN = std::min(iN, i0 + 8);
        // fmt::print("t = {} i0 iN = {} {}\n", time, i0, iN);
    }

    Pos3D result { 0, 0, 0 };
    for (int i = i0; i <= iN; ++i) {
        auto const orbI = vectors[i];

        double weight = 1;
        for (int j = i0; j <= iN; ++j) {
            if (j != i) {
                double const time2 = vectors[j].tp;

                weight *= (time - time2) / (orbI.tp - time2);
            }
        }
        result.x += weight * orbI.xp;
        result.y += weight * orbI.yp;
        result.z += weight * orbI.zp;
    }
    return result;
}
#endif

HD_FUNC
inline double get_zero_doppler_time_newton(double t_guess, Pos3D earth_point, const OrbitPolyInterpolator& orbit)
{
    double t = t_guess;
    constexpr int MAX_ITER = 20;
    constexpr double COND = 1e-8;
    int i = 0;
    for (i = 0; i < MAX_ITER; i++) {
        Pos3D sat_pos = interp_pos(orbit, t);
        Vec3D ep_to_sat = earth_point - sat_pos;
        Vec3D vel_sat = interp_vel(orbit, t);
        double fn = dot(ep_to_sat, vel_sat);
        double dfn = dot(vel_sat, vel_sat);
        // TODO test if acceleration term helps or not in practice
        double dt = fn / dfn;
        t += dt;
        if (std::abs(dt) < COND) {
            break;
        }
    }
    if (i == MAX_ITER) {
        return INVALID_RDR_IDX;
    }
    // fmt::print("iter = {}\n", i);
    return t;
}
