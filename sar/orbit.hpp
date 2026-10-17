// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include "sar_metadata.hpp"

#include "../util/math_utils.hpp"

template <class T, int N>
struct array
{
    // maybe std::array is fine on cuda...
    T arr[N];
};

template <class T, int N>
void cpy(const std::vector<T>& in, array<T, N>& arr)
{
    SARTCPP_ASSERT(in.size() == N);
    std::copy(in.begin(), in.end(), arr.arr);
}

struct XYZ
{
    double x, y, z;
};

struct OrbitPolyInterpolator
{
    static constexpr size_t NP = POLYFIT_DEGREE + 1;
    static constexpr size_t NV = NP - 1;
#if OPI_LAYOUT == 0
    // pos
    array<double, NP> xp;
    array<double, NP> yp;
    array<double, NP> zp;
    // velocity
    array<double, NV> xv;
    array<double, NV> yv;
    array<double, NV> zv;
    // accel
    /*
    array<double, 4> xa;
    array<double, 4> ya;
    array<double, 4> za;
    */
#elif OPI_LAYOUT == 1
    array<XYZ, NP> pos;
    array<XYZ, NV> vel;
#endif
};


#if OPI_LAYOUT == 0
HD_FUNC
inline Pos3D interp_pos(const OrbitPolyInterpolator& orbit, double t)
{
    double xp = polyval(orbit.xp.arr, t);
    double yp = polyval(orbit.yp.arr, t);
    double zp = polyval(orbit.zp.arr, t);
    return {xp, yp, zp};
}

HD_FUNC
inline Vec3D interp_vel(const OrbitPolyInterpolator& orbit, double t)
{
    double xv = polyval(orbit.xv.arr, t);
    double yv = polyval(orbit.yv.arr, t);
    double zv = polyval(orbit.zv.arr, t);
    return {xv, yv, zv};
}
#elif OPI_LAYOUT == 1
HD_FUNC
inline Pos3D interp_pos(const OrbitPolyInterpolator& opi, double t)
{
    double xp = 0.0;
    double yp = 0.0;
    double zp = 0.0;
    const int N = std::size(opi.pos.arr);
    for (int i = 0; i < N; i++)
    {
        xp *= t;
        yp *= t;
        zp *= t;
        xp += opi.pos.arr[i].x;
        yp += opi.pos.arr[i].y;
        zp += opi.pos.arr[i].z;
    }
    return {xp, yp, zp};
}

HD_FUNC
inline Vec3D interp_vel(const OrbitPolyInterpolator& opi, double t)
{
    double xp = 0.0;
    double yp = 0.0;
    double zp = 0.0;
    const int N = std::size(opi.vel.arr);
    for (int i = 0; i < N; i++)
    {
        xp *= t;
        yp *= t;
        zp *= t;
        xp += opi.vel.arr[i].x;
        yp += opi.vel.arr[i].y;
        zp += opi.vel.arr[i].z;
    }
    return {xp, yp, zp};
}
#endif

template <class T, int N>
array<T, N - 1> polyder(array<T, N> in)
{
    array<T, N - 1> out;
    for (int i = 0; i < N - 1; i++)
    {
        out.arr[i] = in.arr[i] * (N - 1 - i);
    }
    return out;
}

template <int N>
array<XYZ, N - 1> polyder(array<XYZ, N> in)
{
    array<XYZ, N - 1> out;
    for (int i = 0; i < N - 1; i++)
    {
        out.arr[i].x = in.arr[i].x * (N - 1 - i);
        out.arr[i].y = in.arr[i].y * (N - 1 - i);
        out.arr[i].z = in.arr[i].z * (N - 1 - i);
    }
    return out;
}

inline void make_poly_interpolator(OrbitPolyInterpolator& opi, const SARMetadata& sar_metadata)
{
    std::vector<double> xp;
    std::vector<double> yp;
    std::vector<double> zp;
    std::vector<double> t;

    //double norm = sar_metadata.osv.back().tp - sar_metadata.osv.front().tp;

    for (const auto& e : sar_metadata.osv)
    {
        xp.push_back(e.xp);
        yp.push_back(e.yp);
        zp.push_back(e.zp);

        t.push_back(e.tp);
    }

    auto poly_xp = polyfit(t, xp, POLYFIT_DEGREE);
    auto poly_yp = polyfit(t, yp, POLYFIT_DEGREE);
    auto poly_zp = polyfit(t, zp, POLYFIT_DEGREE);

#if OPI_LAYOUT == 0
    cpy(poly_xp, opi.xp);
    cpy(poly_yp, opi.yp);
    cpy(poly_zp, opi.zp);

    opi.xv = polyder(opi.xp);
    opi.yv = polyder(opi.yp);
    opi.zv = polyder(opi.zp);
#elif OPI_LAYOUT == 1

    SARTCPP_ASSERT(poly_xp.size() == POLYFIT_DEGREE + 1);
    SARTCPP_ASSERT(poly_yp.size() == poly_xp.size() && poly_zp.size() == poly_xp.size());

    for (int i = 0; i < poly_xp.size(); i++)
    {
        opi.pos.arr[i] = {poly_xp[i], poly_yp[i], poly_zp[i]};
    }

    opi.vel = polyder(opi.pos);
#endif
    // orbit.xa = polyder(poly_xv);
    // orbit.ya = polyder(poly_yv);
    // orbit.za = polyder(poly_zv);
}

inline OSV interpolate_lagrange(const OSV* osv, int n, double tp)
{
    OSV r = {};
    r.tp = tp;
    for (int i = 0; i < n; i++)
    {
        double mult = 1;
        for (int j = 0; j < n; j++)
        {
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
