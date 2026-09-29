// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

// #include "sar_geo.hpp"

#include <array>
#include <fmt/format.h>
#include <gdal/gdal_priv.h>


#include "sar_geo.hpp"
#include "../util/proj_utils.hpp"

struct DEM {
    float* data;
    int x_size;
    int y_size;
    float no_data_value;
    std::array<double, 6> gt;
};

inline void print_spacings_m(const DEM& dem)
{
    double lat_start = dem.gt[3];
    double lon_start = dem.gt[0];
    double pixel_spacing_x = dem.gt[1];
    double pixel_spacing_y = dem.gt[5];

    int mid_y = dem.y_size / 2;
    int mid_x = dem.x_size / 2;

    const double lat = lat_start + mid_y * pixel_spacing_y + 0.5 * pixel_spacing_y;
    const double lon = lon_start + mid_x * pixel_spacing_x + 0.5 * pixel_spacing_x;
    const double lat_p1 = lat + pixel_spacing_y;
    const double lon_p1 = lon + pixel_spacing_x;

    Pos3D p1 = Geo2xyzWgs84(lat, lon, 0);
    Pos3D p2 = Geo2xyzWgs84(lat_p1, lon, 0);
    Pos3D p3 = Geo2xyzWgs84(lat, lon_p1, 0);

    fmt::print("Dem spacing in m wgs84 = {} , {}\n", length(p2 - p1), length(p3 - p1));
}

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

    fmt::print("~~~ DEM ({}) ~~~~\ninfo: sz = ({},{})\ngt = ({},{}) ({},{})\n", path, x_size, y_size, dem.gt[0], dem.gt[3], dem.gt[1], dem.gt[5]);

    int ndv_check = 0;
    double no_data_value = b->GetNoDataValue(&ndv_check);
    if (ndv_check) {
        dem.no_data_value = no_data_value;
        fmt::print("DEM no data value = {}\n", no_data_value);
    }
    auto err = b->RasterIO(GF_Read, 0, 0, x_size, y_size,
        dem.data, x_size, y_size, GDT_Float32,
        0, 0);
    GDALClose(ds);

    SARTCPP_ASSERT(err == CE_None);

    print_spacings_m(dem);
}
