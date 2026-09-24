// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include <fmt/format.h>
#include <gdal/gdal_priv.h>

#include "../util/proj_utils.hpp"

struct DEM {
    float* data;
    int x_size;
    int y_size;
    float no_data_value;
    std::array<double, 6> gt;
};

inline void load_dem(const char* path, DEM& dem)
{
    TimeBlock t("DEM load");
    auto ds = (GDALDataset*)GDALOpen(path, GA_ReadOnly);
    auto b = ds->GetRasterBand(1);
    int x_size = b->GetXSize();
    int y_size = b->GetYSize();
    dem.x_size = x_size;
    dem.y_size = y_size;

    dem.data = new float[x_size * y_size]; // TODO memory leak

    ds->GetGeoTransform(dem.gt.data());

    dem.no_data_value = -12345678.0f;

    fmt::print("DEM({}) info: sz = ({},{}), gt = ({},{}) ({},{})", path, x_size, y_size, dem.gt[0], dem.gt[3], dem.gt[1], dem.gt[5]);

    int ndv_check = 0;
    double no_data_value = b->GetNoDataValue(&ndv_check);
    if (ndv_check) {
        dem.no_data_value = no_data_value;
        fmt::print("DEM no data value = {}\n", no_data_value);
    }
    auto err = b->RasterIO(GF_Read, 0, 0, x_size, y_size,
        dem.data, x_size, y_size, GDT_Float32,
        0, 0);

    SARTCPP_ASSERT(err == CE_None);
}