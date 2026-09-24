// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include "../util/device_raster.hpp"
#include "device_dem.hpp"
#include "sar/sar_metadata.hpp"

void terrain_correct_gpu(const SARMetadata& sar_meta, const DeviceDEM& dem, const DeviceRaster<float>& data_in, DeviceRaster<float>& data_out);
