// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include "s1/s1_metadata.hpp"
#include "sar/sar_metadata.hpp"

#include "../util/device_raster.hpp"

namespace s1 {

void deburst_gpu(const SARMetadata& sar_meta, const S1Metadata& s1_meta, DeviceRaster<float>& in_out);

}