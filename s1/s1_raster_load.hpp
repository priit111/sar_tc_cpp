// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include <optional>

#include "../util/memory_raster.hpp"
#include "s1_metadata.hpp"
#include "../sar/sar_metadata.hpp"

struct SplitParams {
    int cut_start;
    int cut_end;
    SARMetadata* sar_meta;
    s1::S1Metadata* s1_meta;
};

void load_img(const char* path, MemoryRaster<IQ16>& out, std::optional<SplitParams> split_params = std::nullopt);
