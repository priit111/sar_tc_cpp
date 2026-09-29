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
#include "sar/multilook.hpp"

#include "sar/terrain_correction.hpp"
#include "util/arg_parse.hpp"
#include "util/gdal_util.hpp"
#include "util/proj_utils.hpp"

#include "sar/rtc.hpp"

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

    std::optional<SplitParams> split_par = std::nullopt;
    if constexpr (SPLIT_START >=0 && SPLIT_END >= 0) {
        split_par = { SPLIT_START, SPLIT_END, &sar_meta, &s1_meta };
    }
    load_img(sar_meta.raster_path.c_str(), data_in, split_par);

    MemoryRaster<IQ16orF32> calib_arg = data_in.reinterpret_to<IQ16orF32>();
    s1::calibrate(s1_meta, calib_arg);

    MemoryRaster<float> calibrated = calib_arg.reinterpret_to<float>();

    if constexpr (WIF) {
        write_tiff(calibrated, pa.out_root + "cal.tif");
    }

    deburst(sar_meta, s1_meta, calibrated);

    if constexpr (WIF) {
        write_tiff(calibrated, pa.out_root + "deburst.tif");
    }

    if constexpr (EN_RTC) {
        multilook(sar_meta, calibrated);
        if constexpr (WIF) {
            write_tiff(calibrated, pa.out_root + "range_multilook.tif");
        }

        MemoryRaster<float> simulated;
        rtc(sar_meta, dem, simulated, calibrated);

        if constexpr (WIF) {
            write_tiff(simulated, pa.out_root + "sim.tif");
        }
    }

    MemoryRaster<float> tc_out;
    terrain_correct(sar_meta, dem, calibrated, tc_out);

    write_tiff(tc_out, pa.out_path.c_str(), dem.gt, dem.no_data_value);

    return 0;
}
