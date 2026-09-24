// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#include "proj_conf.hpp"

#include <fmt/format.h>
#include <gdal/gdal_priv.h>

#include "s1/s1_calibrate.hpp"
#include "s1/s1_deburst.hpp"
#include "sar/sar_metadata.hpp"

#include "sar/dem.hpp"

#include "s1/s1_parser.hpp"
#include "s1/s1_raster_load.hpp"

#include "sar/sar_geocode.hpp"
#include "util/arg_parse.hpp"
#include "util/gdal_util.hpp"
#include "util/proj_utils.hpp"

int main(int argc, const char* argv[])
{
    fmt::print("sar_tc_cpp version = {}\n", VERSION_STR);
    TimeBlock tb("Total time...");
    if (argc != 6) {
        fmt::print(
            "Missing arguments:\n\n1. S1 SLC directory\n2. DEM\n3. polarisation - "
            "vv or vh\n4. swath - iw1 or iw2 or iw3\n5. output file path\n\n");
        return 1;
    }

    ProgramArgs pa = parse_program_args(argc, argv);

    GDALRegister_GTiff();

    SARMetadata sar_meta = { };
    s1::S1Metadata s1_meta = { };
    s1::parse(pa.s1_dir_path, pa.pol, pa.swath, sar_meta, s1_meta);
    DEM dem = { };
    load_dem(pa.dem_path.c_str(), dem);

    MemoryRaster<IQ16> data_in = { };
    load_img(sar_meta.raster_path.c_str(), data_in);

    MemoryRaster<IQ16orF32> calib_arg = data_in.reinterpret_to<IQ16orF32>();
    s1::calibrate(s1_meta, calib_arg);

    MemoryRaster<float> calibrated = calib_arg.reinterpret_to<float>();

    if constexpr (WIF) {
        write_tiff(calibrated, "/tmp/cal.tiff");
    }

    deburst(sar_meta, s1_meta, calibrated);

    if constexpr (WIF) {
        write_tiff(calibrated, "/tmp/deburst.tiff");
    }

    MemoryRaster<float> tc_out;
    terrain_correct(sar_meta, dem, calibrated, tc_out);

    write_tiff(tc_out, pa.out_path.c_str(), dem.gt);

    return 0;
}
