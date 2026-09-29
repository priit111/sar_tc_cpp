// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include "sar_metadata.hpp"

#include "../util/math_utils.hpp"

template <class T, int N>
struct array { // maybe std::array is fine on cuda...
    T arr[N];
};

template <class T, int N>
void set(const std::vector<T>& in, array<T, N>& arr)
{
    SARTCPP_ASSERT(in.size() == N);
    for (int i = 0; i < N; i++) {
        arr.arr[i] = in[i];
    }
}

struct OrbitPolyInterpolator {
    // pos
    array<double, 6> xp;
    array<double, 6> yp;
    array<double, 6> zp;
    // velocity
    array<double, 5> xv;
    array<double, 5> yv;
    array<double, 5> zv;
    // accel
    /*
    array<double, 4> xa;
    array<double, 4> ya;
    array<double, 4> za;
    */
};

HD_FUNC
inline Pos3D interp_pos(const OrbitPolyInterpolator& orbit, double t)
{
    double xp = polyval(orbit.xp.arr, t);
    double yp = polyval(orbit.yp.arr, t);
    double zp = polyval(orbit.zp.arr, t);
    return { xp, yp, zp };
}

HD_FUNC
inline Vec3D interp_vel(const OrbitPolyInterpolator& orbit, double t)
{
    double xv = polyval(orbit.xv.arr, t);
    double yv = polyval(orbit.yv.arr, t);
    double zv = polyval(orbit.zv.arr, t);
    return { xv, yv, zv };
}

template <class T, int N>
array<T, N - 1> polyder(array<T, N> in)
{
    array<T, N - 1> out;
    for (int i = 0; i < N - 1; i++) {
        out.arr[i] = in.arr[i] * (N - 1 - i);
    }
    return out;
}

inline void make_poly_interpolator(OrbitPolyInterpolator& orbit, const SARMetadata& sar_metadata)
{
    std::vector<double> xp;
    std::vector<double> yp;
    std::vector<double> zp;
    std::vector<double> t;

    for (const auto& e : sar_metadata.osv) {
        xp.push_back(e.xp);
        yp.push_back(e.yp);
        zp.push_back(e.zp);

        t.push_back(e.tp);
    }

    auto poly_xp = polyfit(t, xp, 5);
    auto poly_yp = polyfit(t, yp, 5);
    auto poly_zp = polyfit(t, zp, 5);

    set(poly_xp, orbit.xp);
    set(poly_yp, orbit.yp);
    set(poly_zp, orbit.zp);

    orbit.xv = polyder(orbit.xp);
    orbit.yv = polyder(orbit.yp);
    orbit.zv = polyder(orbit.zp);

    // orbit.xa = polyder(poly_xv);
    // orbit.ya = polyder(poly_yv);
    // orbit.za = polyder(poly_zv);
}

inline OSV interpolate_lagrange(const OSV* osv, int n, double tp)
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