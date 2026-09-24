// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once
#include <thread>

#ifdef __CUDACC__
#define HD_FUNC __host__ __device__
#else
#define HD_FUNC
#endif

inline const char* const VERSION_STR = "0.0.3";

constexpr int SPLIT_START = -1;
constexpr int SPLIT_END = -1;

// performance, algorithm modifications - most will be removed and/or changed to arguments

constexpr int TILE_Y_SIZE = 128; // NB! big memory effects
constexpr int TILE_X_SIZE = 128;

inline int optimal_thread_count()
{
    //return 1;
    return std::thread::hardware_concurrency() / 1;
}

// 0 -> std::mutex
// 1 -> std::atomic
#define TILE_QUEUE_TYPE 0

// 0 -> Bisect
// 1 -> Newton
// TC only, for now RTC is always Newton
#define GEOCODE_MODE 1

#if GEOCODE_MODE == 0

// 0 Lagrange
// 1 Linear
#define TC_POS_INTERP 1
#endif

// GEOCODE_MODE 1
#define ESTIMATE_ZDT 1

// 0 -> NN
// 1 -> bilinear
#define TC_INTERP_MODE 1

// interesting experiment, do restrict pointers matter in this project? Likely not, but keep for now
// #define RESTRICT
#define RESTRICT __restrict__
// #define RESTRICT

// D. Small area distribution
// 0 -> NN
// 1 -> bilinear
#define RTC_INTERP_MODE 1

// 0 -> none
// 1 -> std::atomic_ref<float>
// 2 -> std::atomic_ref<int32_t>
#define RTC_ADD_LOCK 1
