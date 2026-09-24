// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include "../../sar/rtc_types.hpp"

#include "../util/device_raster.hpp"
#include "device_dem.hpp"

void rtc_gpu(const SARMetadata& sar_metadata, RTCParams params, const DeviceDEM& dem, DeviceRaster<float>& out_simulated, DeviceRaster<float>& inout_backscatter);