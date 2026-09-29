// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once


#include <array>
#include <optional>

#include <gdal/gdal_priv.h>

#include "memory_raster.hpp"

inline void write_tiff(const MemoryRaster<float>& in_raster, const char* path,
    std::optional<std::array<double, 6>> gt = std::nullopt, std::optional<float> no_data_value = std::nullopt)
{
    auto ts = TimeStart();

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

inline void write_tiff(const MemoryRaster<float>& in_raster, std::string path,
    std::optional<std::array<double, 6>> gt = std::nullopt, std::optional<float> no_data_value = std::nullopt)
{
    write_tiff(in_raster, path.c_str(), gt, no_data_value);
}
