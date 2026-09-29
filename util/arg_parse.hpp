// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include <string>

#include "proj_utils.hpp"

#include <filesystem>

struct ProgramArgs {
    std::string s1_dir_path;
    std::string dem_path;
    std::string pol;
    std::string swath;
    std::string out_path;
    std::string out_root;
};

inline ProgramArgs parse_program_args(int argc, const char* argv[])
{
    ProgramArgs args;
    SARTCPP_ASSERT(argc == 6);
    args.s1_dir_path = argv[1];
    args.dem_path = argv[2];
    args.pol = argv[3];
    args.swath = argv[4];

    std::filesystem::path p(argv[5]);
    if (std::filesystem::is_directory(p)) {
        p /= "sartcpp_tc.tif";
        args.out_path = p.string();
        args.out_root = p.remove_filename().string();
    }
    else {
        args.out_path = p.string();
        args.out_root = p.remove_filename().string();
    }

    SARTCPP_ASSERT_MSG(std::filesystem::is_directory(args.s1_dir_path), "Input not a dir");
    SARTCPP_ASSERT_MSG(args.pol == "vv" || args.pol == "vh", "invalid pol");
    SARTCPP_ASSERT_MSG(args.swath == "iw1" || args.swath == "iw2" || args.swath == "iw3", "invalid swath");

    fmt::print("ARGS:\nin dir = {}\ndem = {}\npol = {}\nswath = {}\nout_root = {}\nout_path = {}\n",
        args.s1_dir_path, args.dem_path, args.pol, args.swath, args.out_root, args.out_path);

    return args;
}