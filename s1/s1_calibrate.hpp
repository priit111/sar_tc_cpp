// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#pragma once

#include "../util/memory_raster.hpp"
#include "s1_metadata.hpp"

#include <boost/config/no_tr1/memory.hpp>

namespace s1 {
void calibrate(const S1Metadata& s1_meta, MemoryRaster<IQ16orF32>& data_in_out);
}