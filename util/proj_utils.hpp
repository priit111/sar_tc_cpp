#pragma once

// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#include <chrono>
#include <mutex>
#include <optional>
#include <vector>

#include <fmt/format.h>

inline void VerifyFunc(bool condition, const char* filename, const char* function, int line)
{
    if (!condition) {
        fmt::print("SARF_CHECK fail! file = {}  function = {}() line = {}\n", filename, function, line);
        exit(123);
    }
}

#define SARTCPP_ASSERT(a) VerifyFunc(a, __FILE__, __FUNCTION__, __LINE__)

inline auto TimeStart() { return std::chrono::steady_clock::now(); }

inline void TimeStop(std::chrono::steady_clock::time_point start, std::string_view msg)
{
    auto end = std::chrono::steady_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    fmt::print("{} - completed in {} ms\n", msg, diff);
}

struct TimeBlock {
    std::chrono::steady_clock::time_point start;
    std::string msg;
    explicit TimeBlock(std::string msg)
        : msg(std::move(msg))
    {
        start = TimeStart();
    }

    TimeBlock(const TimeBlock&) = delete;
    TimeBlock& operator=(const TimeBlock&) = delete;

    ~TimeBlock()
    {
        TimeStop(start, msg);
    }
};

struct Tile {
    int x_start;
    int x_end;
    int y_start;
    int y_end;
};

inline std::vector<Tile> make_tiles(int target_x_size, int target_y_size, int x_tile_size, int y_tile_size)
{
    std::vector<Tile> tiles;
    int ny = (target_y_size + y_tile_size - 1) / y_tile_size;
    int nx = (target_x_size + x_tile_size - 1) / x_tile_size;
    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            int x_start = x * x_tile_size;
            int x_end = x_start + x_tile_size;
            if (x_end > target_x_size) {
                x_end = target_x_size;
            }
            int y_start = y * y_tile_size;
            int y_end = y_start + y_tile_size;
            if (y_end > target_y_size) {
                y_end = target_y_size;
            }
            tiles.push_back({ x_start, x_end, y_start, y_end });
        }
    }
    return tiles;
}

struct TileQueue {
private:
    std::mutex m_mutex;
    std::vector<Tile> m_vec;
    size_t m_consumed = 0;

public:
    void set_work_blocks(std::vector<Tile>&& work)
    {
        std::lock_guard l(m_mutex);
        m_consumed = 0;
        m_vec = std::move(work);
    }

    std::optional<Tile> get_next_tile()
    {
        std::lock_guard l(m_mutex);
        if (m_consumed == m_vec.size()) {
            return std::nullopt;
        }

        return m_vec[m_consumed++];
    }
};
