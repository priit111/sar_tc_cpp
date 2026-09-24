// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp


#pragma once

#include "s1/s1_metadata.hpp"

#include "../util/device_raster.hpp"

namespace s1 {
void calibrate_gpu(const S1Metadata& s1_meta, DeviceRaster<IQ16orF32>& data_in_out);
}