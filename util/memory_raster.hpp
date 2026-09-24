// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include <cstdint>
#include <type_traits>

#include "fmt/format.h"

#include "proj_utils.hpp"

struct IQ16 {
    int16_t i;
    int16_t q;
};

/*
template <>
struct fmt::formatter<IQ16> : fmt::formatter<std::string> {
    auto format(const IQ16& iq16, format_context& ctx) const {
        return fmt::formatter<std::string>::format(
            fmt::format("({}, {})", iq16.i, iq16.q), ctx);
    }
};
*/

union IQ16orF32 {
    float f32;
    IQ16 iq16;
};

static_assert(sizeof(IQ16orF32) == 4);

template <class T>
struct MemoryRaster {
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, IQ16> || std::is_same_v<T, IQ16orF32>);
    T* m_data;
    int m_x_size;
    int m_y_size;

    void init(int x_size, int y_size)
    {
        SARTCPP_ASSERT(m_data == nullptr);
        m_data = new T[x_size * y_size];
        m_x_size = x_size;
        m_y_size = y_size;
    }

    size_t total_size() const
    {
        return static_cast<size_t>(m_x_size) * m_y_size;
    }

    T get_pixel(int x, int y) const
    {
        return m_data[y * m_x_size + x];
    }

    void resize_y(int new_y)
    {
        SARTCPP_ASSERT(new_y < m_y_size);
        m_y_size = new_y;
    }

    MemoryRaster()
        : m_data(nullptr)
        , m_x_size(0)
        , m_y_size { 0 }
    {
    }

    MemoryRaster(T* ext, int x_size, int y_size)
        : m_data(ext)
        , m_x_size(x_size)
        , m_y_size(y_size)

    {
    }
    MemoryRaster(const MemoryRaster&) = delete;
    MemoryRaster& operator=(const MemoryRaster&) = delete;
    MemoryRaster(MemoryRaster&& other) noexcept
    {
        m_data = other.m_data;
        m_x_size = other.m_x_size;
        m_y_size = other.m_y_size;
        other.m_data = nullptr;
        other.m_x_size = 0;
        other.m_y_size = 0;
    }
    MemoryRaster& operator=(MemoryRaster&& other) noexcept
    {
        m_data = other.m_data;
        m_x_size = other.m_x_size;
        m_y_size = other.m_y_size;
        other.m_data = nullptr;
        other.m_x_size = 0;
        other.m_y_size = 0;
        return *this;
    }
    ~MemoryRaster()
    {
        delete[] m_data;
    }

    template <class U>
    MemoryRaster<U> reinterpret_to()
    {
        static_assert(sizeof(U) == sizeof(T));
        MemoryRaster<U> other;
        other.m_data = reinterpret_cast<U*>(m_data); // TODO use/read about start_lifetime_as
        other.m_x_size = m_x_size;
        other.m_y_size = m_y_size;

        m_x_size = 0;
        m_y_size = 0;
        m_data = nullptr;
        return other;
    }
};
