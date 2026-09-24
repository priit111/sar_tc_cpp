// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include <string>

#include "proj_utils.hpp"

struct ProgramArgs {
    std::string s1_dir_path;
    std::string dem_path;
    std::string pol;
    std::string swath;
    std::string out_path;
};

inline ProgramArgs parse_program_args(int argc, const char* argv[])
{
    ProgramArgs args;
    SARTCPP_ASSERT(argc == 6);
    args.s1_dir_path = argv[1];
    args.dem_path = argv[2];
    args.pol = argv[3];
    args.swath = argv[4];
    args.out_path = argv[5];

    SARTCPP_ASSERT(args.pol == "vv" || args.pol == "vh");
    SARTCPP_ASSERT(args.swath == "iw1" || args.swath == "iw2" || args.swath == "iw3");

    return args;
}