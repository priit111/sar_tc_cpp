// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#include "rtc.hpp"

#include "orbit.hpp"

#include "../proj_conf.hpp"

#include "sar_geo.hpp"

// 0 -> NN
// 1 -> bilinear
#define RTC_INTERP_MODE 1

// 0 -> none
// 1 -> atomic_ref

#define RTC_ADD_LOCK 0

struct Dimensions {
    int dem_x_size;
    int dem_y_size;

    int az_max_idx;
    int rg_max_idx;
};

struct Args {
    double az_spacing;
    double rg_spacing;

    double line_time_interval;
    double slant_range_first_sample;
    double lon_start;
    double lat_start;
    double pixel_spacing_x;
    double pixel_spacing_y;
    float no_data_value;
    OrbitPolyInterpolator opi;
};

float blerp_dem(const float* dem, double x, double y, int x_size, float no_data_value)
{
    int x_i = static_cast<int>(x);
    int y_i = static_cast<int>(y);
    float dx = x - x_i;
    float dy = y - y_i;
    float alt00 = dem[y_i * x_size + x_i];
    float alt10 = dem[y_i * x_size + x_i + 1];
    float alt01 = dem[(y_i + 1) * x_size + x_i];
    float alt11 = dem[(y_i + 1) * x_size + x_i + 1];
    if (alt00 == no_data_value || alt01 == no_data_value || alt10 == no_data_value || alt11 == no_data_value) {
        return no_data_value;
    }
    return blerp(alt00, alt10, alt01, alt11, dx, dy);
}

std::mutex s_mutex;

void fetch_add(float* ptr, float val)
{
#if RTC_ADD_LOCK == 0
    *ptr += val;
#endif
#if RTC_ADD_LOCK == 1
    std::atomic_ref<float> ar(*ptr);
    ar.fetch_add(val, std::memory_order_relaxed); // NB not deterministic due to fp rounding?
#endif
}

void calc_rtc(TileQueue* tq, const float* dem, Dimensions io_dims, Args args, float* sim_out)
{
    const int dem_y_size = io_dims.dem_y_size;
    const int dem_x_size = io_dims.dem_x_size;

    const int sim_x_size = io_dims.rg_max_idx;
    const int sim_y_size = io_dims.az_max_idx;
    const int max_az_idx = io_dims.az_max_idx;
    const int max_rg_idx = io_dims.rg_max_idx;
    const double rg_m_to_pix = 1.0 / args.rg_spacing;

    const double dem_x_step = 1.0 / RTC_X_OVERSAMP;
    const double dem_y_step = 1.0 / RTC_Y_OVERSAMP;

    while (1) {
        double t_guess = 5;
        std::optional<Tile> next_tile = tq->get_next_tile();
        if (!next_tile.has_value()) {
            break;
        }
        Tile tile = next_tile.value();
        double y = tile.y_start;
        const double y_end = tile.y_end;
        const double x_end = tile.x_end;
        while (y < y_end) {

            double x = tile.x_start;
            while (x < x_end) {
                // fmt::print("x = {} , y = {}\n", x,y);
                const double lat = args.lat_start + y * args.pixel_spacing_y + 0.5 * args.pixel_spacing_y;
                const double lon = args.lon_start + x * args.pixel_spacing_x + 0.5 * args.pixel_spacing_x;

                float altitude = blerp_dem(dem, x, y, dem_x_size, args.no_data_value);

                if (altitude == args.no_data_value) {
                    x += dem_x_step;
                    continue;
                }
                float alt_E00 = altitude;
                float alt_E10 = blerp_dem(dem, x + dem_x_step, y, dem_x_size, args.no_data_value); // dem[y * dem_x_size + x + 1];
                float alt_E01 = blerp_dem(dem, x, y + dem_y_step, dem_x_size, args.no_data_value); // dem[(y + 1) * dem_x_size + x];
                float alt_E11 = blerp_dem(dem, x + dem_x_step, y + dem_y_step, dem_x_size, args.no_data_value); // dem[(y + 1) * dem_x_size + x + 1];
                if (alt_E10 == args.no_data_value || alt_E01 == args.no_data_value || alt_E11 == args.no_data_value) {
                    x += dem_x_step;
                    continue;
                }

                const double lat_Ex0 = lat;
                const double lat_Ex1 = lat + dem_y_step * args.pixel_spacing_y;
                const double lon_E0x = lon;
                const double lon_E1x = lon + dem_x_step * args.pixel_spacing_x;

                auto E00 = Geo2xyzWgs84(lat_Ex0, lon_E0x, alt_E00);
                auto E10 = Geo2xyzWgs84(lat_Ex0, lon_E1x, alt_E10);

                auto E01 = Geo2xyzWgs84(lat_Ex1, lon_E0x, alt_E01);
                auto E11 = Geo2xyzWgs84(lat_Ex1, lon_E1x, alt_E11);

                double az_time = get_zero_doppler_time_newton(t_guess, E00, args.opi);
                t_guess = az_time;
                double az_idx = az_time / args.line_time_interval;
                if (az_idx < 0 || az_idx > io_dims.az_max_idx) {
                    x += dem_x_step;

                    continue;
                }

                Pos3D sat_pos = interp_pos(args.opi, az_time);

                double slant_range = distance(E00, sat_pos);
                double rg_idx = (slant_range - args.slant_range_first_sample) * rg_m_to_pix;
                if (rg_idx < 0 || rg_idx >= max_rg_idx) {
                    x += dem_x_step;

                    continue;
                }

                // computeIlluminatedArea
                // TerrainFlatteningOp.java#L1314
                // for now the D.Small geometry calcs taken from SNAP
                Vec3D slant_dir = sat_pos - E00;
                slant_dir = normalize(slant_dir);

                double t00s = dot(from_origin(E00), slant_dir);
                double t01s = dot(from_origin(E01), slant_dir);

                double t10s = dot(from_origin(E10), slant_dir);
                double t11s = dot(from_origin(E11), slant_dir);

                Pos3D p00 = { E00.x - t00s * slant_dir.x, E00.y - t00s * slant_dir.y, E00.z - t00s * slant_dir.z };
                Pos3D p01 = { E01.x - t01s * slant_dir.x, E01.y - t01s * slant_dir.y, E01.z - t01s * slant_dir.z };
                Pos3D p10 = { E10.x - t10s * slant_dir.x, E10.y - t10s * slant_dir.y, E10.z - t10s * slant_dir.z };
                Pos3D p11 = { E11.x - t11s * slant_dir.x, E11.y - t11s * slant_dir.y, E11.z - t11s * slant_dir.z };

                // compute distances between projected points
                double p00p01 = distance(p00, p01);
                double p00p10 = distance(p00, p10);
                double p11p01 = distance(p11, p01);
                double p11p10 = distance(p11, p10);
                double p10p01 = distance(p10, p01);

                double h1 = 0.5 * (p00p01 + p00p10 + p10p01);
                double h2 = 0.5 * (p11p01 + p11p10 + p10p01);

                double A = sqrt(h1 * (h1 - p00p01) * (h1 - p00p10) * (h1 - p10p01)) + sqrt(h2 * (h2 - p11p01) * (h2 - p11p10) * (h2 - p10p01));

                int rg_idx_i = std::round(rg_idx);
                int az_idx_i = std::round(az_idx);
#if RTC_INTERP_MODE == 0
                int rg_idx_i = std::round(rg_idx);
                int az_idx_i = std::round(az_idx);

                float* ptr00 = &sim_out[az_idx_i * sim_x_size + rg_idx_i];
                fetch_add(ptr00, A);
#elif RTC_INTERP_MODE == 1
                float Wr = rg_idx - rg_idx_i;
                float Wa = az_idx - az_idx_i;
                float Wcr = 1 - Wr;
                float Wca = 1 - Wa;

                float A00 = Wcr * Wca * A;
                float A10 = Wr * Wca * A;
                float A01 = Wcr * Wa * A;
                float A11 = Wr * Wa * A;

                float* ptr00 = &sim_out[az_idx_i * sim_x_size + rg_idx_i];
                float* ptr01 = &sim_out[(az_idx_i + 1) * sim_x_size + rg_idx_i];
                float* ptr10 = &sim_out[az_idx_i * sim_x_size + rg_idx_i + 1];
                float* ptr11 = &sim_out[(az_idx_i + 1) * sim_x_size + rg_idx_i + 1];
                fetch_add(ptr00, A00);
                fetch_add(ptr10, A10);
                fetch_add(ptr01, A01);
                fetch_add(ptr11, A11);

                x += dem_x_step;

#endif
            }
            y += dem_y_step;
        }
    }
}

#if 0
void calc_rtc2(TileQueue* tq, const float* dem, Dimensions io_dims, Args args, float* sim_out)
{
    const int dem_y_size = io_dims.dem_y_size;
    const int dem_x_size = io_dims.dem_x_size;

    const int sim_x_size = io_dims.rg_max_idx;
    const int sim_y_size = io_dims.az_max_idx;
    const int max_az_idx = io_dims.az_max_idx;
    const int max_rg_idx = io_dims.rg_max_idx;
    const double rg_m_to_pix = 1.0 / args.rg_spacing;


    while (1) {
        double t_guess = 5;
        std::optional<Tile> next_tile = tq->get_next_tile();
        if (!next_tile.has_value()) {
            break;
        }
        Tile tile = next_tile.value();
        double y = tile.y_start;
        const double y_end = tile.y_end;
        const double x_end = tile.x_end;

        for (int y = tile.y_start; y < y_end; y++) {

            for (int x = tile.x_start; x < x_end; x++) {
                // fmt::print("x = {} , y = {}\n", x,y);
                const double lat = args.lat_start + y * args.pixel_spacing_y + 0.5 * args.pixel_spacing_y;
                const double lon = args.lon_start + x * args.pixel_spacing_x + 0.5 * args.pixel_spacing_x;

                float altitude = dem[y * dem_x_size + x];

                if (altitude == args.no_data_value) {
                    continue;
                }
                float alt_E00 = altitude;
                float alt_E10 = dem[y * dem_x_size + x + 1];
                float alt_E01 = dem[(y + 1) * dem_x_size + x];
                float alt_E11 = dem[(y + 1) * dem_x_size + x + 1];
                if (alt_E10 == args.no_data_value || alt_E01 == args.no_data_value || alt_E11 == args.no_data_value) {
                    continue;
                }

                const double lat_Ex0 = lat;
                const double lat_Ex1 = lat +  args.pixel_spacing_y;
                const double lon_E0x = lon;
                const double lon_E1x = lon +  args.pixel_spacing_x;

                auto E00 = Geo2xyzWgs84(lat_Ex0, lon_E0x, alt_E00);
                auto E10 = Geo2xyzWgs84(lat_Ex0, lon_E1x, alt_E10);

                auto E01 = Geo2xyzWgs84(lat_Ex1, lon_E0x, alt_E01);
                auto E11 = Geo2xyzWgs84(lat_Ex1, lon_E1x, alt_E11);

                Pos3D cent_pos = center(E00, E10, E11);

                Vec3D e1 = E10 - E00;
                Vec3D e2 = E11 - E00;


                Vec3D n = cross(e1, e2);

                if (dot(n, from_origin(E00)) < 0.0) {
                    n = n * -1.0;
                }

                double az_time = get_zero_doppler_time_newton(t_guess, cent_pos, args.opi);

                Pos3D sat_pos = interp_pos(args.opi, az_time);
                Vec3D look = sat_pos - cent_pos;
                look = normalize(look);
                auto l_hat = look;
                    double area = 0.5 * dot(n, l_hat);
                if (area < 0.0) {
                    area = 0.0;
                }
                double az_t_E00 = get_zero_doppler_time_newton(t_guess, E00, args.opi);
                double az_t_E10 = get_zero_doppler_time_newton(t_guess, E10, args.opi);
                double az_t_E01 = get_zero_doppler_time_newton(t_guess, E01, args.opi);
                double az_t_E11 = get_zero_doppler_time_newton(t_guess, E11, args.opi);

                double az_time = get_zero_doppler_time_newton(t_guess, E00, args.opi);
                t_guess = az_time;
                double az_idx = az_time / args.line_time_interval;
                if (az_idx < 0 || az_idx > io_dims.az_max_idx) {
                    continue;
                }

                Pos3D sat_pos = interp_pos(args.opi, az_time);

                double slant_range = distance(E00, sat_pos);
                double rg_idx = (slant_range - args.slant_range_first_sample) * rg_m_to_pix;
                if (rg_idx < 0 || rg_idx >= max_rg_idx) {

                    continue;
                }

                // computeIlluminatedArea
                // TerrainFlatteningOp.java#L1314
                // for now the D.Small papers geometry calcs taken from SNAP
                Vec3D slant_dir = sat_pos - E00;
                slant_dir = normalize(slant_dir);

                double t00s = dot(from_origin(E00), slant_dir);
                double t01s = dot(from_origin(E01), slant_dir);

                double t10s = dot(from_origin(E10), slant_dir);
                double t11s = dot(from_origin(E11), slant_dir);

                Pos3D p00 = { E00.x - t00s * slant_dir.x, E00.y - t00s * slant_dir.y, E00.z - t00s * slant_dir.z };
                Pos3D p01 = { E01.x - t01s * slant_dir.x, E01.y - t01s * slant_dir.y, E01.z - t01s * slant_dir.z };
                Pos3D p10 = { E10.x - t10s * slant_dir.x, E10.y - t10s * slant_dir.y, E10.z - t10s * slant_dir.z };
                Pos3D p11 = { E11.x - t11s * slant_dir.x, E11.y - t11s * slant_dir.y, E11.z - t11s * slant_dir.z };

                // compute distances between projected points
                double p00p01 = distance(p00, p01);
                double p00p10 = distance(p00, p10);
                double p11p01 = distance(p11, p01);
                double p11p10 = distance(p11, p10);
                double p10p01 = distance(p10, p01);

                double h1 = 0.5 * (p00p01 + p00p10 + p10p01);
                double h2 = 0.5 * (p11p01 + p11p10 + p10p01);

                double A = sqrt(h1 * (h1 - p00p01) * (h1 - p00p10) * (h1 - p10p01)) + sqrt(h2 * (h2 - p11p01) * (h2 - p11p10) * (h2 - p10p01));

                int rg_idx_i = std::round(rg_idx);
                int az_idx_i = std::round(az_idx);
#if RTC_INTERP_MODE == 0
                int rg_idx_i = std::round(rg_idx);
                int az_idx_i = std::round(az_idx);

                float* ptr00 = &sim_out[az_idx_i * sim_x_size + rg_idx_i];
                fetch_add(ptr00, A);
#elif RTC_INTERP_MODE == 1
                float Wr = rg_idx - rg_idx_i;
                float Wa = az_idx - az_idx_i;
                float Wcr = 1 - Wr;
                float Wca = 1 - Wa;

                float A00 = Wcr * Wca * A;
                float A10 = Wr * Wca * A;
                float A01 = Wcr * Wa * A;
                float A11 = Wr * Wa * A;

                float* ptr00 = &sim_out[az_idx_i * sim_x_size + rg_idx_i];
                float* ptr01 = &sim_out[(az_idx_i + 1) * sim_x_size + rg_idx_i];
                float* ptr10 = &sim_out[az_idx_i * sim_x_size + rg_idx_i + 1];
                float* ptr11 = &sim_out[(az_idx_i + 1) * sim_x_size + rg_idx_i + 1];
                fetch_add(ptr00, A00);
                fetch_add(ptr10, A10);
                fetch_add(ptr01, A01);
                fetch_add(ptr11, A11);

#endif
            }
        }
    }
}
#endif

struct RTCFinalizeArgs {
    int sim_x_size;
    int sim_y_size;
    float incidence_begin;
    float incidence_end;
    float ref_area;
};

void finalize_rtc(RTCFinalizeArgs args, float* simulated, float* backscatter_beta)
{
    int y_size = args.sim_y_size;
    int x_size = args.sim_x_size;
    constexpr float D2R = M_PI / 180.0f;

    // for now just assume linear incidence angle change in range
    // close enough and almost no impact on calculation time
    const float k_incidence = (args.incidence_end - args.incidence_begin) / x_size;
    for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
            float& A = simulated[y * x_size + x];
            float inci = args.incidence_begin + k_incidence * x;
            float a = args.ref_area / tan(D2R * inci);
            A /= a;
            if (backscatter_beta) {
                if (A == 0) {
                    backscatter_beta[y * x_size + x] = 0;
                } else {
                    backscatter_beta[y * x_size + x] /= A;
                }
            }
        }
    }
}

void rtc(const SARMetadata& sar_metadata, const DEM& dem, MemoryRaster<float>& out_simulated, MemoryRaster<float>& inout_backscatter)
{
    fmt::print("Dsmall start\n");
    TimeBlock tb("RTC D.Small");
    out_simulated.init(sar_metadata.range_size, sar_metadata.azimuth_size);
    out_simulated.clear();

    Dimensions io_dims = { };
    io_dims.dem_x_size = dem.x_size;
    io_dims.dem_y_size = dem.y_size;
    io_dims.az_max_idx = sar_metadata.azimuth_size;
    io_dims.rg_max_idx = sar_metadata.range_size;

    Args args = { };
    args.line_time_interval = sar_metadata.line_time_interval;

    make_poly_interpolator(args.opi, sar_metadata);

    args.slant_range_first_sample = sar_metadata.slant_range_first_sample;
    args.lat_start = dem.gt[3];
    args.lon_start = dem.gt[0];
    args.pixel_spacing_x = dem.gt[1];
    args.pixel_spacing_y = dem.gt[5];
    args.no_data_value = dem.no_data_value;

    args.rg_spacing = sar_metadata.range_spacing;
    args.az_spacing = sar_metadata.azimuth_spacing;

    const int n_threads = optimal_thread_count();
    std::vector<std::thread> thread_vec;

    constexpr size_t y_block_sz = TILE_Y_SIZE;
    constexpr size_t x_block_sz = TILE_X_SIZE;

    fmt::print("RTC conf N THREADS = {} tile = [{} , {}], interp = {}\n",
        n_threads, x_block_sz, y_block_sz, RTC_INTERP_MODE);

    if (n_threads == 1) {
        TileQueue tq;
        tq.set_work_blocks({ { 0, io_dims.dem_x_size, 0, io_dims.dem_y_size } });
        calc_rtc(&tq, dem.data, io_dims, args, out_simulated.m_data);
    } else {
        TileQueue tq;
        tq.set_work_blocks(make_tiles(io_dims.dem_x_size, io_dims.dem_y_size, x_block_sz, y_block_sz));

        for (int i = 0; i < n_threads; i++) {

            std::thread t(calc_rtc, &tq, dem.data, io_dims, args, out_simulated.m_data);
            thread_vec.push_back(std::move(t));
        }
        for (auto& t : thread_vec) {
            t.join();
        }
    }

    {
        TimeBlock tb2("finalize rtc");
        RTCFinalizeArgs fargs = { };
        fargs.ref_area = sar_metadata.azimuth_spacing * sar_metadata.range_spacing;
        fargs.incidence_begin = sar_metadata.incidence_angle_begin;
        fargs.incidence_end = sar_metadata.incidence_angle_end;

        fargs.sim_x_size = sar_metadata.range_size;
        fargs.sim_y_size = sar_metadata.azimuth_size;
        finalize_rtc(fargs, out_simulated.m_data, inout_backscatter.m_data);
    }
}