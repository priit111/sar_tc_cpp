// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#include "proj_conf.hpp"

#include <thread>

#include <fmt/format.h>
#include <gdal/gdal_priv.h>

#include "s1/s1_calibrate_gpu.hpp"
#include "s1/s1_deburst_gpu.hpp"
#include "sar/sar_metadata.hpp"

#include "sar/dem.hpp"

#include "s1/s1_parser.hpp"
#include "s1/s1_raster_load.hpp"

#include "sar/sar_geocode_gpu.hpp"
#include "util/arg_parse.hpp"
#include "util/gdal_util.hpp"
#include "util/proj_utils.hpp"

#include "sar/device_dem.hpp"

#include "util/cuda_util.hpp"
#include "util/device_raster.h"

int main(int argc, const char* argv[])
{
    fmt::print("sar_tc_cuda version = {}\n", VERSION_STR);

    TimeBlock tb("Total time...");
    if (argc != 6) {
        fmt::print(
            "Missing arguments:\n\n1. S1 SLC directory\n2. DEM\n3. polarisation - "
            "vv or vh\n4. swath - iw1 or iw2 or iw3\n5. output file path\n\n");
        return 1;
    }

    ProgramArgs pa = parse_program_args(argc, argv);

    std::thread cuda_init([]() { cudaFree(nullptr); });
    cuda_init.detach();

    GDALRegister_GTiff();

    SARMetadata sar_meta = { };
    s1::S1Metadata s1_meta = { };
    s1::parse(pa.s1_dir_path, pa.pol, pa.swath, sar_meta, s1_meta);
    DEM dem = { };
    load_dem(pa.dem_path.c_str(), dem);

    MemoryRaster<IQ16> h_data_in = { };
    load_img(sar_meta.raster_path.c_str(), h_data_in);

    DeviceRaster<IQ16> d_data_in = { };

    d_data_in.init(h_data_in.m_x_size, h_data_in.m_y_size);

    h2d_cpy(d_data_in.m_d_data, h_data_in.m_data, h_data_in.total_size());

    DeviceRaster<IQ16orF32> calib_arg = d_data_in.reinterpret_to<IQ16orF32>();

    s1::calibrate_gpu(s1_meta, calib_arg);
    DeviceRaster<float> calibrated = calib_arg.reinterpret_to<float>();
    if (WIF) {
        auto hr = d2h_raster(calibrated);
        write_tiff(hr, "/tmp/s1_calib_cuda.tif");
    }

    s1::deburst_gpu(sar_meta, s1_meta, calibrated);

    if (WIF) {
        auto hr = d2h_raster(calibrated);
        write_tiff(hr, "/tmp/deburst_cuda.tif");
    }
    DeviceDEM d_dem;
    d_dem.init(dem.data, dem.x_size, dem.y_size, dem.no_data_value, dem.gt);
    DeviceRaster<float> d_tc_out;
    terrain_correct_gpu(sar_meta, d_dem, calibrated, d_tc_out);

    //* reusing host DEM memory...
    d2h_cpy(dem.data, d_tc_out.m_d_data, d_tc_out.total_size());
    MemoryRaster<float> dem_steal(dem.data, dem.x_size, dem.y_size);
    write_tiff(dem_steal, pa.out_path.c_str(), dem.gt);
    return 0;
}
