// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include "../util/device_raster.hpp"
#include "../../sar/sar_metadata.hpp"

void multilook_gpu(SARMetadata& sar_meta, DeviceRaster<float>& memory_raster, int range_multilook);