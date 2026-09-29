// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include "../proj_conf.hpp"

#include <cmath>
#include <type_traits>
#include <vector>

/*
q12    q22
 ______
 |    |
 |____|
q11   q21
 */

template <class T>
HD_FUNC
inline T blerp(
    T q11, T q21,
    T q12, T q22,
    T x, T y)
{
    // x and y are assumed to be in [0, 1]
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>);
    T bot = q11 + (q21 - q11) * x;
    T top = q12 + (q22 - q12) * x;

    return bot + (top - bot) * y;
}

// constexpr float TEST = blerp(11.0f, 22.0f, 33.0f, 44.0f, 0.6f, 0.8f);
// static_assert(TEST == 35.2f); // FMA rounding... TODO


/*

            | i    j    k  |
  a x b  =  | ax   ay   az |
            | bx   by   bz |
 */

struct Vec3D {
    double x, y, z;
};
struct Pos3D {
    double x, y, z;
};

HD_FUNC
inline Vec3D operator-(Pos3D a, Pos3D b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

HD_FUNC
inline Vec3D operator* (Vec3D a, double k)
{
    a.x *= k;
    a.y *= k;
    a.z *= k;
    return a;
}

HD_FUNC
inline Vec3D from_origin(Pos3D pos)
{
    return  pos - Pos3D{0, 0, 0};
}

HD_FUNC
inline double length(Vec3D a)
{
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

HD_FUNC
inline Vec3D normalize(Vec3D a)
{
    double len = length(a);
    return {a.x / len, a.y / len, a.z / len};
}

HD_FUNC
inline double distance(Pos3D a, Pos3D b)
{
    return length(a-b);
}

HD_FUNC
inline Vec3D cross(Vec3D a, Vec3D b)
{
    return {(a.y * b.z - a.z * b.y), -(a.x * b.z - b.x * a.z),(a.x * b.y - b.x * a.y)};
}

HD_FUNC
inline double dot(Vec3D a, Vec3D b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

HD_FUNC
inline Pos3D center(Pos3D a, Pos3D b, Pos3D c)
{
    constexpr double mul = 1/3.0;
    return {(a.x + b.x + c.x) * mul, (a.y + b.y + c.y) * mul, (a.z + b.z + c.z) * mul };
}


template<int N>
HD_FUNC
double polyval(const double(&p)[N], double x)
{
    double val = 0.0;
    for (int i = 0; i < N; i++) {
        val *= x;
        val += p[i];
    }
    return val;
}

std::vector<double> polyfit(const std::vector<double>& x, const std::vector<double>& y, int degree);

