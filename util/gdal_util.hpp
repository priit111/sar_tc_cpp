// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include <gdal/gdal_priv.h>
#include <optional>

#include "memory_raster.hpp"

inline void write_tiff(const MemoryRaster<float>& in_raster, const char* path, std::optional<std::array<double, 6>> gt = std::nullopt)
{
    auto ts = TimeStart();

    GDALRegister_GTiff();
    int h = in_raster.m_y_size;
    int w = in_raster.m_x_size;
    auto ds = GetGDALDriverManager()->GetDriverByName("gtiff")->Create(path, w, h, 1, GDT_Float32, nullptr);

    SARTCPP_ASSERT(ds);
    auto b = ds->GetRasterBand(1);
    auto err = b->RasterIO(GF_Write, 0, 0, w, h, in_raster.m_data, w, h, GDT_Float32, 0, 0);
    SARTCPP_ASSERT(err == CE_None);

    if (gt.has_value()) {
        ds->SetGeoTransform(gt->data());
    }
    GDALClose(ds);
    TimeStop(ts, fmt::format("File @ {}", path));
}

inline void write_direct_tiff(const MemoryRaster<float>& in_raster, const char* path, std::optional<std::array<double, 6>> gt = std::nullopt)
{
    auto ts = TimeStart();
    int h = in_raster.m_y_size;
    int w = in_raster.m_x_size;
    auto ds = GetGDALDriverManager()->GetDriverByName("gtiff")->Create(path, w, h, 1, GDT_Float32, nullptr);

    SARTCPP_ASSERT(ds);
    auto b = ds->GetRasterBand(1);
    int x_block_sz = 0;
    int y_block_sz = 0;
    b->GetBlockSize(&x_block_sz, &y_block_sz);
    SARTCPP_ASSERT(x_block_sz == w);
    SARTCPP_ASSERT(y_block_sz == 1);

    float* data_in = const_cast<float*>(in_raster.m_data);
    for (int y = 0; y < h; y++) {
        // float*
        auto err = b->WriteBlock(0, y, data_in);
        SARTCPP_ASSERT(err == CE_None);
        data_in += x_block_sz;
    }

    if (gt.has_value()) {
        ds->SetGeoTransform(gt->data());
    }
    GDALClose(ds);

    TimeStop(ts, fmt::format("Striped file @ {}", path));
}