// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include "../util/memory_raster.hpp"
#include "sar_metadata.hpp"

void multilook(SARMetadata& sar_meta, MemoryRaster<float>& memory_raster, int range_multilook);
