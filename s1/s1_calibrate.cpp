// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#include "s1_calibrate.hpp"

#include <thread>

namespace {
void beta_calibrate(TileQueue* tq, IQ16orF32* data, int x_size, float beta_coeff)
{

    while (1) {
        std::optional<Tile> next_tile = tq->get_next_tile();
        if (!next_tile.has_value()) {
            break;
        }
        Tile tile = next_tile.value();
        int y_start = tile.y_start;
        int y_end = tile.y_end;
        int x_start = tile.x_start;
        int x_end = tile.x_end;
        for (int y = y_start; y < y_end; y++) {
            for (int x = x_start; x < x_end; x++) {
                size_t idx = y * x_size + x;
                IQ16 iq16 = data[idx].iq16;
                float I = iq16.i;
                float Q = iq16.q;
                float cal = (I * I + Q * Q) * beta_coeff;
                data[idx].f32 = cal;
            }
        }
    }
}
} // namespace

namespace s1 {
void calibrate(const S1Metadata& s1_meta, MemoryRaster<IQ16orF32>& data_in_out)
{
    float beta = s1_meta.calib.front().beta.front(); // an hack for now, assume beta = constant

    const int n_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> thread_vec;
    int y_size = data_in_out.m_y_size;
    int x_size = data_in_out.m_x_size;
    int y_step = (y_size / n_threads) + 1;

    TileQueue tq;
    tq.set_work_blocks(make_tiles(x_size, y_size, x_size, 100));
    for (int i = 0; i < n_threads; i++) {
        std::thread t(beta_calibrate, &tq, data_in_out.m_data, x_size, 1.0f / (beta * beta));

        thread_vec.push_back(std::move(t));
    }

    for (auto& t : thread_vec) {
        t.join();
    }
}
} // namespace s1