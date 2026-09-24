#pragma once


#include <cuda_runtime.h>

#include "util/memory_raster.hpp"
#include "cuda_util.hpp"

static_assert(sizeof(IQ16orF32) == 4);

template <class T>
struct DeviceRaster {
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, IQ16> || std::is_same_v<T, IQ16orF32>);
    T* m_d_data;
    int m_x_size;
    int m_y_size;

    void init(int x_size, int y_size)
    {
        SARTCPP_ASSERT(m_d_data == nullptr);
        T* ptr;

        auto err = cudaMalloc(&ptr, static_cast<size_t>(x_size) * y_size * sizeof(T));
        CUDA_VERIFY_CHECK(err);
        m_d_data = ptr;
        m_x_size = x_size;
        m_y_size = y_size;
    }

    size_t total_size() const
    {
        return static_cast<size_t>(m_x_size) * m_y_size;
    }

    T get_pixel(int x, int y) const
    {
        return m_d_data[y * m_x_size + x];
    }

    void resize_y(int new_y)
    {
        SARTCPP_ASSERT(new_y < m_y_size);
        m_y_size = new_y;
    }

    DeviceRaster()
        : m_d_data(nullptr)
        , m_x_size(0)
        , m_y_size { 0 }
    {
    }
    DeviceRaster(const DeviceRaster&) = delete;
    DeviceRaster& operator=(const DeviceRaster&) = delete;
    DeviceRaster(DeviceRaster&& other) noexcept
    {
        m_d_data = other.m_d_data;
        m_x_size = other.m_x_size;
        m_y_size = other.m_y_size;
        other.m_d_data = nullptr;
        other.m_x_size = 0;
        other.m_y_size = 0;
    }
    DeviceRaster& operator=(DeviceRaster&& other) noexcept
    {
        m_d_data = other.m_d_data;
        m_x_size = other.m_x_size;
        m_y_size = other.m_y_size;
        other.m_d_data = nullptr;
        other.m_x_size = 0;
        other.m_y_size = 0;
        return *this;
    }
    ~DeviceRaster()
    {
        cudaFree(m_d_data);
    }

    template <class U>
    DeviceRaster<U> reinterpret_to()
    {
        static_assert(sizeof(U) == sizeof(T));
        DeviceRaster<U> other;
        other.m_d_data = reinterpret_cast<U*>(m_d_data); // TODO use/read about start_lifetime_as
        other.m_x_size = m_x_size;
        other.m_y_size = m_y_size;

        m_x_size = 0;
        m_y_size = 0;

        m_d_data = nullptr;
        return other;
    }
};


inline MemoryRaster<float> d2h_raster(const DeviceRaster<float>& d_raster)
{
    MemoryRaster<float> h;
    h.init(d_raster.m_x_size, d_raster.m_y_size);
    d2h_cpy(h.m_data, d_raster.m_d_data, h.total_size());
    return h;
}