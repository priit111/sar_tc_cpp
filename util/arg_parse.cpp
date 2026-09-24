// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#include "arg_parse.hpp"


ProgramArgs parse_program_args(int argc, const char* argv[])
{
    ProgramArgs args;
    SARTCPP_ASSERT(argc >= 6);
    args.s1_dir_path = argv[1];
    args.dem_path = argv[2];
    args.pol = argv[3];
    args.swath = argv[4];

    std::filesystem::path p(argv[5]);
    if (std::filesystem::is_directory(p)) {
        p /= "sartcpp_tc.tif";
        args.out_path = p.string();
        args.out_root = p.remove_filename().string();
    } else {
        args.out_path = p.string();
        args.out_root = p.remove_filename().string();
    }

    args.rtc_en = false;
    args.rtc_x_upsamp = args.rtc_y_upsamp = 0;
    args.range_multilook = 1;
    args.wif = false;

    if (argc >= 6) {
        for (int i = 6; i < argc; i++) {
            fmt::print("OPT ARG = {}\n", argv[i]);
            std::string_view sv = argv[i];
            if (sv.starts_with("--rtc=")) {
                float x = -1;
                float y = -1;
                if (sscanf(sv.data(), "--rtc=%f,%f", &x, &y) == 2) {
                    args.rtc_x_upsamp = x;
                    args.rtc_y_upsamp = y;
                }

                SARTCPP_ASSERT_MSG(x >= 1.0 && y >= 1.0, "Arg parse error");
                args.rtc_en = true;
            } else if (sv.starts_with("--wif")) {
                args.wif = true;
            } else if (sv.starts_with("--ml")) {
                int r = -1;
                if (sscanf(sv.data(), "--ml=%d", &r) == 1) {
                    args.range_multilook = r;
                }
                SARTCPP_ASSERT_MSG(r > 1, "Arg parse_error");
            } else {
                SARTCPP_ASSERT_MSG(false, "Arg parse error");
            }
        }
    }

    SARTCPP_ASSERT_MSG(std::filesystem::is_directory(args.s1_dir_path), "Input not a dir");
    SARTCPP_ASSERT_MSG(args.pol == "vv" || args.pol == "vh", "invalid pol");
    SARTCPP_ASSERT_MSG(args.swath == "iw1" || args.swath == "iw2" || args.swath == "iw3", "invalid swath");

    fmt::print("ARGS:\nin dir = {}\ndem = {}\npol = {}\nswath = {}\nout_root = {}\nout_path = {}\n",
        args.s1_dir_path, args.dem_path, args.pol, args.swath, args.out_root, args.out_path);
    fmt::print("RTC en = {}\n", args.rtc_en);
    if (args.rtc_en) {
        fmt::print("RTC upsamp = {} , {}\nrange multilook = {}\n", args.rtc_x_upsamp, args.rtc_y_upsamp, args.range_multilook);
    }

    return args;
}