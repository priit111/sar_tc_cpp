// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#include "s1_raster_load.hpp"

#include <fmt/format.h>
#include <gdal/gdal_priv.h>

#include "../util/memory_raster.hpp"

void load_img(const char* path, MemoryRaster<IQ16>& out)
{
    TimeBlock tb(fmt::format("raster = {} load", path));
    GDALDataset* ds = (GDALDataset*)GDALOpen(path, GA_ReadOnly);

    SARTCPP_ASSERT(ds);
    auto b = ds->GetRasterBand(1);

    int x_size = b->GetXSize();
    int y_size = b->GetYSize();
    out.init(x_size, y_size);
    auto err = b->RasterIO(GF_Read, 0, 0, x_size, y_size,
        out.m_data, x_size, y_size, GDT_CInt16,
        0, 0);

    // TODO readblock api usage

    SARTCPP_ASSERT(err == CE_None);

    GDALClose(ds);
}
