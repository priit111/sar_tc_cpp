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
    bool rtc_en;
    double rtc_x_upsamp;
    double rtc_y_upsamp;
    int range_multilook;
    bool wif;
};

ProgramArgs parse_program_args(int argc, const char* argv[]);
