// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#include "s1_raster_load.hpp"

#include <fmt/format.h>
#include <gdal/gdal_priv.h>

#include "../util/memory_raster.hpp"
#include "../util/proj_utils.hpp"

void load_img(const char* path, MemoryRaster<IQ16>& out, std::optional<SplitParams> split_params)
{
    TimeBlock tb(fmt::format("raster = {} load", path));
    GDALDataset* ds = (GDALDataset*)GDALOpen(path, GA_ReadOnly);

    SARTCPP_ASSERT(ds);
    auto b = ds->GetRasterBand(1);

    int x_size = b->GetXSize();
    int y_size = b->GetYSize();

    int y_start = 0;
    int y_end = y_size;
    if (split_params.has_value()) {
        SARMetadata& sar_meta = *split_params->sar_meta;
        s1::S1Metadata& s1_meta = *split_params->s1_meta;
        int cut_start = split_params->cut_start;
        int cut_end = split_params->cut_end;
        int lpb = s1_meta.lines_per_burst;
        int n_bursts = s1_meta.bursts.size();
        SARTCPP_ASSERT(lpb * n_bursts == y_size);
        SARTCPP_ASSERT(cut_end > cut_start);
        SARTCPP_ASSERT(cut_start >= 0 && cut_end < n_bursts);
        y_start = cut_start * lpb;
        y_end = cut_end * lpb;
        y_size = y_end - y_start;
        fmt::print("Using bursts [{} , {})\n", cut_start, cut_end);
        fmt::print("y offsets = {} {}\n", y_start, y_end);

        s1_meta.bursts.resize(cut_end);
        s1_meta.bursts.erase(s1_meta.bursts.begin(), s1_meta.bursts.begin() + cut_start);

        fmt::print("Burst after = {}\n", s1_meta.bursts.size());

        sar_meta.cut_y(y_size, s1_meta.bursts.front().az_time);
    }

    out.init(x_size, y_size);
    auto err = b->RasterIO(GF_Read, 0, y_start, x_size, y_size,
        out.m_data, x_size, y_size, GDT_CInt16,
        0, 0);

    // TODO readblock api usage ?

    SARTCPP_ASSERT(err == CE_None);

    GDALClose(ds);
}
