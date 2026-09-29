// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once
#include <thread>


#ifdef __CUDACC__
#define HD_FUNC __host__ __device__
#else
#define HD_FUNC
#endif



inline const char* const VERSION_STR = "0.0.2";

constexpr bool WIF = false;
constexpr int SPLIT_START = -1;
constexpr int SPLIT_END = -1;



// performance, algorithm modifications - most will be removed and/or changed to arguments

constexpr int TILE_Y_SIZE = 1; // NB! big memory effects
constexpr int TILE_X_SIZE = 1024;

inline int optimal_thread_count()
{
    return std::thread::hardware_concurrency() / 2;
}




// 0 -> Bisect
// 1 -> Newton
#define GEOCODE_MODE 1

#if GEOCODE_MODE == 0

// 0 Lagrange
// 1 Linear
#define TC_POS_INTERP 1
#endif

// 0 -> NN
// 1 -> bilinear
#define TC_INTERP_MODE 1


// interesting experiment, do restrict pointers matter in this project? Likely not, but keep for now
// #define RESTRICT
#define RESTRICT __restrict__


// WIP, very first attempt at D. Small RTC, CPU only for now
constexpr bool EN_RTC = false;
constexpr int RTC_RANGE_MULTILOOK = 4;
// D. Small area distribution
// 0 -> NN
// 1 -> bilinear
#define RTC_INTERP_MODE 1

constexpr double RTC_Y_OVERSAMP = 1.5;
constexpr double RTC_X_OVERSAMP = 1.5;