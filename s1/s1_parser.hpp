// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once
#include <string>

#include "../sar/sar_metadata.hpp"
#include "s1_metadata.hpp"

namespace s1 {
bool parse(std::string s1_dir_path, std::string pol, std::string swath, SARMetadata& sar_meta, S1Metadata& s1_meta);
}
