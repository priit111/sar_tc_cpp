// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp


#pragma once

#include "sar_metadata.hpp"
#include "../util/memory_raster.hpp"

void multilook(SARMetadata& sar_meta, MemoryRaster<float>& memory_raster);
