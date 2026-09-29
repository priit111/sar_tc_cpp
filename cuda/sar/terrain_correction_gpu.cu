
#include "terrain_correction_gpu.hpp"

#include "../util/cuda_util.hpp"


#include "sar/sar_metadata.hpp"

#include "sar/sar_geo.hpp"


// #include "../util/math_utils.hpp"

// ~~~~~~~~~~~
// NB very copy paste from sar_geocode.cpp
// TODO better host-device code sharing
// ~~~~~~~~~~~


constexpr bool MASK_NO_DATA_VALUE = true;

// CONF END

namespace {


__global__ void RunTC(const float* RESTRICT data_in, const float* RESTRICT dem, TCDimensions io_dims, TCArgs args, float* RESTRICT data_out)
{
    const int x = threadIdx.x + blockDim.x * blockIdx.x;
    const int y = threadIdx.y + blockDim.y * blockIdx.y;

    if (x > io_dims.out_x_size || y > io_dims.out_y_size) {
        return;
    }

    // std::vector<GeoPos3D> pos_vec(args.sat_pos_az, args.sat_pos_az + io_dims.in_y_size);
    // std::vector<Velocity3D> vel_vec(args.sat_vel_az, args.sat_vel_az + io_dims.in_y_size);
    // args.sat_pos_az = pos_vec.data();
    // args.sat_vel_az = vel_vec.data();
    const int out_x_size = io_dims.out_x_size;
    const int out_y_size = io_dims.out_y_size;
    const int in_x_size = io_dims.in_x_size;
    const int in_y_size = io_dims.in_y_size;
    const double rg_m_to_pix = 1 / args.range_spacing;
#if GEOCODE_MODE == 1
    double guess_t = in_y_size * args.line_time_interval / 2; // TODO prediction sharing between kernels, easy on CPU, needs thinking on cuda
#endif

    const double lat = args.lat_start + y * args.pixel_spacing_y + 0.5 * args.pixel_spacing_y;
    const double lon = args.lon_start + x * args.pixel_spacing_x + 0.5 * args.pixel_spacing_x;

    const int out_idx = x + y * out_x_size;

    float altitude = dem[y * out_x_size + x];
    if (altitude == args.no_data_value) {
        if (MASK_NO_DATA_VALUE) {
            data_out[out_idx] = 0;
        } else {
            altitude = 0;
        }
    }

    auto earth_point = Geo2xyzWgs84(lat, lon, altitude);

#if GEOCODE_MODE == 0
    const double az_time = GetEarthPointZeroDopplerTime(args.line_time_interval, args.wavelength, earth_point,
        io_dims.in_y_size, args.sat_pos_az, args.sat_vel_az);
#elif GEOCODE_MODE == 1
    double az_time = get_zero_doppler_time_newton(guess_t, earth_point, args.opi);
    double az_interp_time = az_time;
    guess_t = az_time;
    // az_time2 *= args.line_time_interval;

    if (x == out_x_size / 2 && y == out_y_size / 2) {
        // fmt::print("AZ times = {} {} {}\n", az_time, az_time2, az_interp_time);
    }

#endif
    // const int out_idx = x + y * out_x_size;
    if (az_time == INVALID_RDR_IDX) {
        data_out[out_idx] = 0;
        return;
    }
    // fmt::print("VALID !\n");
    const double az_idx = az_time / args.line_time_interval;

    if (az_idx < 0.0 || az_idx >= in_y_size - 1) {
        data_out[out_idx] = 0;
        return;
    }

    // slant range
#if GEOCODE_MODE == 0
#if TC_POS_INTERP == 0
    // SNAP interpolate OSV for each point
    GeoPos3D sat_pos = GetPosition(az_time, args.osv, args.n_osv);
#elif TC_POS_INTERP == 1
    // linearly interpolate between precalculated positions? Faster with almost no precision loss?
    // TODO investigate
    double interp_dt = az_idx - std::round(az_idx);
    int az_idx_i = static_cast<int>(az_idx);
    auto pos0 = args.sat_pos_az[az_idx_i];
    auto pos1 = args.sat_pos_az[az_idx_i + 1];
    Pos3D sat_pos = { };
    sat_pos.x = pos0.x + interp_dt * (pos1.x - pos0.x);
    sat_pos.y = pos0.y + interp_dt * (pos1.y - pos0.y);
    sat_pos.z = pos0.z + interp_dt * (pos1.z - pos0.z);
#endif
#elif GEOCODE_MODE == 1
    Pos3D sat_pos = interp_pos(args.opi, az_interp_time);

#endif
    // range idx
    double slant_range = distance(earth_point, sat_pos);
    double rg_idx = (slant_range - args.slant_range_first_sample) * rg_m_to_pix;
    if (rg_idx < 0 || rg_idx >= in_x_size) {
        data_out[out_idx] = 0;
        return;
    }

#if TC_INTERP_MODE == 0
    int rg_idx_nn = std::round(rg_idx);
    int az_idx_nn = std::round(az_idx);

    size_t in_idx = rg_idx_nn + in_x_size * az_idx_nn;
    data_out[out_idx] = data_in[in_idx];
#elif TC_INTERP_MODE == 1
    int rg_idx_int = static_cast<int>(rg_idx);
    int az_idx_int = static_cast<int>(az_idx);
    size_t in_idx = rg_idx_int + in_x_size * az_idx_int;
    // TODO corner idx check
    float q11 = data_in[in_idx]; // x floor, y floor
    float q21 = data_in[in_idx + 1]; // x ceil, y floor
    float q12 = data_in[in_idx + in_x_size]; // x floor, y ceil
    float q22 = data_in[in_idx + in_x_size + 1]; // x ceil, y ceil
    // float
    float mx = rg_idx - rg_idx_int;
    float my = az_idx - az_idx_int;
    float val = blerp(q11, q21, q12, q22, mx, my);
    data_out[out_idx] = val;
#else
#error "..."
#endif
}
} // namespace

void terrain_correct_gpu(const SARMetadata& sar_meta, const DeviceDEM& dem, const DeviceRaster<float>& data_in, DeviceRaster<float>& data_out)
{
    TimeBlock t("TC gpu");
    int range_size = dem.x_size;
    int azimuth_size = dem.y_size;

    data_out.init(range_size, azimuth_size);

    TCArgs args = { };
    args.lat_start = dem.gt[3];
    args.lon_start = dem.gt[0];
    args.pixel_spacing_x = dem.gt[1];
    args.pixel_spacing_y = dem.gt[5];

    args.no_data_value = dem.no_data_value;

    args.line_time_interval = sar_meta.line_time_interval;
    args.wavelength = sar_meta.wavelength;
    args.range_spacing = sar_meta.range_spacing;
    args.slant_range_first_sample = sar_meta.slant_range_first_sample;

#if GEOCODE_MODE == 0
    std::vector<Pos3D> pos;
    std::vector<Vec3D> vel;
    for (int i = 0; i < data_in.m_y_size; i++) {
        auto osv = InterpolateOrbit(sar_meta.osv.data(), sar_meta.osv.size(), sar_meta.calc_az_tp(i));
        pos.push_back({ osv.xp, osv.yp, osv.zp });
        vel.push_back({ osv.xv, osv.yv, osv.zv });
    }

    auto osv = sar_meta.osv;

    osv.clear();

    double last_line_time = sar_meta.calc_az_tp(azimuth_size - 1);
    for (double t = -5.0; t < last_line_time + 5.0; t += 1.0) {
        auto iosv = InterpolateOrbit(sar_meta.osv.data(), sar_meta.osv.size(), t);
        osv.push_back(iosv);
    }

    DeviceBuffer<Pos3D> d_pos;
    d_pos.init(pos.size());
    d_pos.cpy_from_host(pos.data());

    DeviceBuffer<Vec3D> d_vel;
    d_vel.init(vel.size());
    d_vel.cpy_from_host(vel.data());

    DeviceBuffer<OSV> d_osv;

    d_osv.init(osv.size());
    d_osv.cpy_from_host(osv.data());

    args.sat_pos_az = d_pos.data();
    args.sat_vel_az = d_vel.data();
    args.osv = d_osv.data();
    args.n_osv = d_osv.size();
#elif GEOCODE_MODE == 1

    make_poly_interpolator(args.opi, sar_meta);
#endif


    TCDimensions io_dims = { };
    io_dims.in_x_size = data_in.m_x_size;
    io_dims.in_y_size = data_in.m_y_size;
    io_dims.out_x_size = range_size;
    io_dims.out_y_size = azimuth_size;

    dim3 block_size(16, 16);
    dim3 grid_size((io_dims.out_x_size + 15) / 16, (io_dims.out_y_size + 15) / 16);

    RunTC<<<grid_size, block_size>>>(data_in.m_d_data, dem.d_data.m_d_data, io_dims, args, data_out.m_d_data);

    CUDA_VERIFY_CHECK(cudaDeviceSynchronize()); // not needed at this point, but helps with debugging
    CUDA_VERIFY_CHECK(cudaGetLastError());
}
