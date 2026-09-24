// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#include "rtc_gpu.hpp"

#include "../util/device_raster.hpp"

#include "../../sar/orbit.hpp"

#include "../../proj_conf.hpp"

#include "../../sar/sar_geo.hpp"

namespace {
__device__ float blerp_dem(const float* dem, double x, double y, int x_size, float no_data_value)
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

__device__ void fetch_add(float* ptr, float val)
{
#if RTC_ADD_LOCK == 0
    *ptr += val;
#elif RTC_ADD_LOCK == 1
    atomicAdd(ptr, val);
#elif RTC_ADD_LOCK == 2
    int32_t* p = reinterpret_cast<int32_t*>(ptr);
    int32_t val_i = val * 1000;
    atomicAdd(p, val_i);
#else
#error "..."
#endif
}

struct CudaExtraArgs {
    double dem_x_step;
    double dem_y_step;
    ZDT_Estimate zdt_estimate;
};


__global__ void calc_rtc_kernel(const float* RESTRICT dem, RTCDimensions io_dims, RTCCalcArgs args, float* RESTRICT sim_out, CudaExtraArgs extra_args)
{



    const int dem_y_size = io_dims.dem_y_size;
    const int dem_x_size = io_dims.dem_x_size;

    const int sim_x_size = io_dims.rg_max_idx;
    const int sim_y_size = io_dims.az_max_idx;
    const int max_az_idx = io_dims.az_max_idx;
    const int max_rg_idx = io_dims.rg_max_idx;


    const int x_start = threadIdx.x + blockDim.x * blockIdx.x;
    const int y_start = threadIdx.y + blockDim.y * blockIdx.y;

    if (x_start >= io_dims.dem_x_size - 2 || y_start >= io_dims.dem_y_size - 2) {
        return;
    }


    const double rg_m_to_pix = args.rg_m_to_pix;

    //const double dem_x_step = 1.0 / args.x_oversamp_ratio;
    //const double dem_y_step = 1.0 / args.y_oversamp_ratio;
    double dem_x_step = extra_args.dem_x_step;
    double dem_y_step = extra_args.dem_y_step;



    double t_guess = 5;
#if ESTIMATE_ZDT == 1
    t_guess = estimate_zdt(extra_args.zdt_estimate, x_start, y_start);
#endif
    double y = y_start;
    const double y_end = y_start + 1;;
    const double x_end = x_start + 1;
    while (y < y_end) {
        double x = x_start;
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
            float alt_E10 = blerp_dem(dem, x + dem_x_step, y, dem_x_size, args.no_data_value);
            // dem[y * dem_x_size + x + 1];
            float alt_E01 = blerp_dem(dem, x, y + dem_y_step, dem_x_size, args.no_data_value);
            // dem[(y + 1) * dem_x_size + x];
            float alt_E11 = blerp_dem(dem, x + dem_x_step, y + dem_y_step,
                dem_x_size, args.no_data_value);
            // dem[(y + 1) * dem_x_size + x + 1];
            if (alt_E00 == args.no_data_value || alt_E10 == args.no_data_value || alt_E01 == args.no_data_value || alt_E11 == args.no_data_value) {
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
            double az_idx = az_time * args.az_t_to_pix;
            if (az_idx < 0 || az_idx > max_az_idx) {
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

            Vec3D slant_dir = sat_pos - E00;
            slant_dir = normalize(slant_dir);
#if 1
            Vec3D n1 = cross((E01 - E00), (E10 - E00));
            Vec3D n2 = cross((E01 - E11), (E10 - E11));

            // fmt::print("dot n1 n2 = {}\n", dot(n1, n2));

            if (dot(from_origin(E00), n1) < 0) {
                n1 = n1 * -1.0;
            }

            if (dot(from_origin(E11), n2) < 0) {
                n2 = n2 * -1.0;
            }

            //TODO sign handling...
            double A1 = 0.5 * (dot(n1, slant_dir));
            double A2 = 0.5 * (dot(n2, slant_dir));

            double A = A1 + A2;

#endif
#if 0
            // TerrainFlatteningOp.java#L1314 for comparision
            //
            double t00s = dot(from_origin(E00), slant_dir);
            double t01s = dot(from_origin(E01), slant_dir);

            double t10s = dot(from_origin(E10), slant_dir);
            double t11s = dot(from_origin(E11), slant_dir);

            Pos3D p00 = {E00.x - t00s * slant_dir.x, E00.y - t00s * slant_dir.y, E00.z - t00s * slant_dir.z};
            Pos3D p01 = {E01.x - t01s * slant_dir.x, E01.y - t01s * slant_dir.y, E01.z - t01s * slant_dir.z};
            Pos3D p10 = {E10.x - t10s * slant_dir.x, E10.y - t10s * slant_dir.y, E10.z - t10s * slant_dir.z};
            Pos3D p11 = {E11.x - t11s * slant_dir.x, E11.y - t11s * slant_dir.y, E11.z - t11s * slant_dir.z};

            // compute distances between projected points
            double p00p01 = distance(p00, p01);
            double p00p10 = distance(p00, p10);
            double p11p01 = distance(p11, p01);
            double p11p10 = distance(p11, p10);
            double p10p01 = distance(p10, p01);

            double h1 = 0.5 * (p00p01 + p00p10 + p10p01);
            double h2 = 0.5 * (p11p01 + p11p10 + p10p01);

            double A = sqrt(h1 * (h1 - p00p01) * (h1 - p00p10) * (h1 - p10p01)) + sqrt(
                h2 * (h2 - p11p01) * (h2 - p11p10) * (h2 - p10p01));


            double Asnap1 = sqrt(h1 * (h1 - p00p01) * (h1 - p00p10) * (h1 - p10p01));
            double Asnap2 = sqrt(h2 * (h2 - p11p01) * (h2 - p11p10) * (h2 - p10p01));
            //fmt::print("{} {} , {} {}\n", Asnap1, Asnap2, A1, A2);
#endif
#if RTC_INTERP_MODE == 0
            int rg_idx_i = std::round(rg_idx);
            int az_idx_i = std::round(az_idx);

            float* ptr00 = &sim_out[az_idx_i * sim_x_size + rg_idx_i];
            fetch_add(ptr00, A);
#elif RTC_INTERP_MODE == 1
            int rg_idx_i = std::floor(rg_idx);
            int az_idx_i = std::floor(az_idx);
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
            //printf("KERNEL ADD %f %f\n", x, y);
#endif
            x += dem_x_step;

        }
        y += dem_y_step;
    }

        //printf("DONE! %d %d\n", x_start, y_start);
}

struct RTCFinalizeArgs {
    int sim_x_size;
    int sim_y_size;
    float incidence_begin;
    float incidence_end;
    float ref_area;
};

__global__ void finalize_rtc_kernel(RTCFinalizeArgs args, float* RESTRICT simulated,
    float* RESTRICT backscatter_beta)
{

    const int x = threadIdx.x + blockDim.x * blockIdx.x;
    const int y = threadIdx.y + blockDim.y * blockIdx.y;

    if (x >= args.sim_x_size || y >= args.sim_y_size) {
        return;
    }
    int y_size = args.sim_y_size;
    int x_size = args.sim_x_size;
    constexpr float D2R = M_PI / 180.0f;

    // for now just assume linear incidence angle change in range
    // close enough and minimal impact on calculation time
    const float k_incidence = (args.incidence_end - args.incidence_begin) / x_size;

    float& A = simulated[y * x_size + x];
#if RTC_ADD_LOCK == 2
    int32_t val;
    memcpy(&val, &A, 4);
    float t = val * 0.001;
    A = t;
#endif
    //float inci = args.incidence_begin + k_incidence * x;
    float a = args.ref_area; // / tan(D2R * inci);

    A /= a;
    if (backscatter_beta) {
        if (A < 0.015) {
            backscatter_beta[y * x_size + x] = 0;
        } else {
            backscatter_beta[y * x_size + x] /= A;
        }
    }

}
}


void rtc_gpu(const SARMetadata& sar_metadata, RTCParams params, const DeviceDEM& dem, DeviceRaster<float>& out_simulated, DeviceRaster<float>& inout_backscatter)
{
    fmt::print("D.Small RTC gpu start\n");
    fmt::print("X,Y oversamp ratio = {} , {}, interp mode = {}\n", params.x_oversamp_ratio,
        params.y_oversamp_ratio, RTC_INTERP_MODE);
    TimeBlock tb("RTC D.Small");
    out_simulated.init(sar_metadata.range_size, sar_metadata.azimuth_size);
    out_simulated.zero_fill();

    RTCDimensions io_dims = {};
    io_dims.dem_x_size = dem.x_size;
    io_dims.dem_y_size = dem.y_size;
    io_dims.az_max_idx = sar_metadata.azimuth_size;
    io_dims.rg_max_idx = sar_metadata.range_size;

    RTCCalcArgs args = {};
    args.az_t_to_pix = 1.0 / sar_metadata.line_time_interval;

    make_poly_interpolator(args.opi, sar_metadata);

    args.slant_range_first_sample = sar_metadata.slant_range_first_sample;
    args.lat_start = dem.gt[3];
    args.lon_start = dem.gt[0];
    args.pixel_spacing_x = dem.gt[1];
    args.pixel_spacing_y = dem.gt[5];
    args.x_oversamp_ratio = params.x_oversamp_ratio;
    args.y_oversamp_ratio = params.y_oversamp_ratio;
    args.no_data_value = dem.no_data_value;

    args.rg_m_to_pix = 1.0 / sar_metadata.range_spacing;

    // blerp_dem_test(dem, args, io_dims);

    CudaExtraArgs extra = {};

    extra.dem_x_step = 1.0 / params.x_oversamp_ratio;
    extra.dem_y_step = 1.0 / params.y_oversamp_ratio;

    extra.zdt_estimate = create_zdt_estimate(args.opi, dem.gt, dem.x_size, dem.y_size);

    {
        dim3 block_size(16, 16);
        dim3 grid_size((io_dims.dem_x_size + 15) / 16, (io_dims.dem_y_size + 15) / 16);

        calc_rtc_kernel<<<grid_size, block_size>>>(dem.d_data.m_d_data, io_dims, args, out_simulated.m_d_data, extra);

        CUDA_VERIFY_CHECK(cudaDeviceSynchronize()); // not needed at this point, but helps with debugging
        CUDA_VERIFY_CHECK(cudaGetLastError());
    }
    
    {
        TimeBlock tb2("finalize rtc");
        RTCFinalizeArgs fargs = {};
        fargs.ref_area = sar_metadata.azimuth_spacing * sar_metadata.range_spacing;
        fargs.incidence_begin = sar_metadata.incidence_angle_begin;
        fargs.incidence_end = sar_metadata.incidence_angle_end;

        fargs.sim_x_size = sar_metadata.range_size;
        fargs.sim_y_size = sar_metadata.azimuth_size;

        dim3 block_size(16, 16);
        dim3 grid_size((fargs.sim_x_size + 15) / 16, (fargs.sim_y_size + 15) / 16);

        finalize_rtc_kernel<<<grid_size, block_size>>>(fargs, out_simulated.m_d_data, inout_backscatter.m_d_data);
        CUDA_VERIFY_CHECK(cudaDeviceSynchronize()); // not needed at this point, but helps with debugging
        CUDA_VERIFY_CHECK(cudaGetLastError());
        //finalize_rtc(fargs, out_simulated.m_data, inout_backscatter.m_data);
    }
}
