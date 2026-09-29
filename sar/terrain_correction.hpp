// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include "../util/memory_raster.hpp"
#include "dem.hpp"
#include "sar_metadata.hpp"

void terrain_correct(const SARMetadata& sar_meta, const DEM& dem, const MemoryRaster<float>& data_in, MemoryRaster<float>& data_out);