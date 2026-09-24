// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

/*
q12    q22
 ______
 |    |
 |____|
q11   q21
 */
template <class T>
constexpr T blerp(
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