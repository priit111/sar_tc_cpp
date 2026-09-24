#pragma once

#include <stdlib.h>

#include <cuda_runtime.h>

#include "util/proj_utils.hpp"

inline void VerifyCudaFunc(cudaError_t e, const char* filename, const char* function, int line)
{
    if (e != cudaSuccess) {
        fmt::print("CUDA_CHECK fail! e = {} file = {}  function = {}() line = {}\n", (int)e, filename, function, line);
        exit(123);
    }
}

#define CUDA_VERIFY_CHECK(a) VerifyCudaFunc(a, __FILE__, __FUNCTION__, __LINE__)

template <class T>
void h2d_cpy(T* d_dest, const T* h_src, size_t n_elem)
{
    static_assert(std::is_pod_v<T>, "POD only on gpu!");
    auto err = cudaMemcpy(d_dest, h_src, n_elem * sizeof(T), cudaMemcpyHostToDevice);
    CUDA_VERIFY_CHECK(err);
}

template <class T>
void d2h_cpy(T* h_dst, const T* d_src, size_t n_elem)
{
    static_assert(std::is_pod_v<T>, "POD only on gpu!");
    auto err = cudaMemcpy(h_dst, d_src, n_elem * sizeof(T), cudaMemcpyDeviceToHost);
    CUDA_VERIFY_CHECK(err);
}

template <class T>
struct DeviceBuffer {
    size_t m_n_elem;
    T* m_d_data;
    DeviceBuffer()

        : m_n_elem(0)
        , m_d_data(nullptr)
    {
    }
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    ~DeviceBuffer()
    {
        cudaFree(m_d_data);
    }

    T* data()
    {
        return m_d_data;
    }

    const T* data() const
    {
        return m_d_data;
    }

    void init(size_t n_elem)
    {
        T* p = nullptr;
        cudaMalloc(&p, n_elem * sizeof(T));

        m_n_elem = n_elem;
        m_d_data = p;
    }

    size_t size() const
    {
        return m_n_elem;
    }

    void cpy_from_host(const T* h_src)
    {
        h2d_cpy(m_d_data, h_src, m_n_elem);
    }

    void cpy_to_host(const T* h_dst)
    {
        d2h_cpy(h_dst, m_d_data, m_n_elem);
    }
};