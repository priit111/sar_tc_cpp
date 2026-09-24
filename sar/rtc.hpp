// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include "../util/memory_raster.hpp"
#include "dem.hpp"
#include "sar_metadata.hpp"
#include "orbit.hpp"
#include "rtc_types.hpp"

void rtc(const SARMetadata& sar_metadata, RTCParams params, const DEM& dem, MemoryRaster<float>& out_simulated, MemoryRaster<float>& inout_backscatter);