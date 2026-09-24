
#include "sar_geocode_gpu.hpp"

#include "../util/cuda_util.hpp"

#include <thread>

#include "sar/sar_metadata.hpp"

#include <boost/range/algorithm/equal.hpp>

// #include "../util/math_utils.hpp"

// ~~~~~~~~~~~
// NB very copy paste from sar_geocode.cpp
// TODO better host-device code sharing
// ~~~~~~~~~~~

// WIP TC implementation testing configuration
// CONF START

// 0 -> NN
// 1 -> bilinear
#define INTERP_MODE 1

// #define RESTRICT
#define RESTRICT __restrict__

// 0 Lagrange
// 1 Linear
#define TC_POS_INTERP 1

constexpr bool ZERO_DOPPLER_INTERCEPT = false;

constexpr bool MASK_NO_DATA_VALUE = true;

// CONF END

namespace {
struct GeoPos3D {
    double x;
    double y;
    double z;
};
struct Velocity3D {
    double x;
    double y;
    double z;
};

namespace WGS84 { // NOLINT
    constexpr double A = 6378137.0; // m
    constexpr double B = 6356752.3142451794975639665996337; // 6356752.31424518; // m
    constexpr double FLAT_EARTH_COEF = 1.0 / ((A - B) / A); // 298.257223563;
    constexpr double E2 = 2.0 / FLAT_EARTH_COEF - 1.0 / (FLAT_EARTH_COEF * FLAT_EARTH_COEF);
    constexpr double EP2 = E2 / (1 - E2);
} // namespace WGS84

constexpr double DTOR = M_PI / 180.0;
constexpr double RTOD = 180 / M_PI;

__device__ double Distance(GeoPos3D xyz1, GeoPos3D xyz2)
{
    double dx = xyz1.x - xyz2.x;
    double dy = xyz1.y - xyz2.y;
    double dz = xyz1.z - xyz2.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

__device__ GeoPos3D Geo2xyzWgs84(double latitude, double longitude, double altitude)
{
    double const lat = latitude * DTOR;
    double const lon = longitude * DTOR;

    double sin_lat;
    double cos_lat;
    sincos(lat, &sin_lat, &cos_lat);

    double const sinLat = sin_lat;

    double const N = (WGS84::A / sqrt(1.0 - WGS84::E2 * sinLat * sinLat));
    double const NcosLat = (N + altitude) * cos_lat;

    double sin_lon;
    double cos_lon;
    sincos(lon, &sin_lon, &cos_lon);

    GeoPos3D r = { };
    r.x = NcosLat * cos_lon;
    r.y = NcosLat * sin_lon;
    r.z = (N + altitude - WGS84::E2 * N) * sinLat;
    return r;
}

__device__ double CalcDopplerFrequency(GeoPos3D earth_point, GeoPos3D sensor_pos, Velocity3D sensor_vel, double wavelength)
{
    const auto dx = earth_point.x - sensor_pos.x;
    const auto dy = earth_point.y - sensor_pos.y;
    const auto dz = earth_point.z - sensor_pos.z;
    const auto distance = sqrt(dx * dx + dy * dy + dz * dz);

    return 2.0 * (sensor_vel.x * dx + sensor_vel.y * dy + sensor_vel.z * dz) / (distance * wavelength);
}

struct Dimensions {
    int in_x_size;
    int in_y_size;
    int out_x_size;
    int out_y_size;
};
struct Args {
    // geo lat/lon
    double lon_start;
    double lat_start;
    double pixel_spacing_y;
    double pixel_spacing_x;

    float no_data_value;

    // geometry of the acquisition
    double line_time_interval;
    double range_spacing;
    double slant_range_first_sample;

    // radar
    double wavelength;

    // OSV
    const OSV* osv;
    int n_osv;

    // LUTs for the zero Doppler binary search
    const GeoPos3D* sat_pos_az;
    const Velocity3D* sat_vel_az;
};

constexpr double INVALID_AZ_TIME = -1234567e7;
__device__ double GetEarthPointZeroDopplerTime(double line_time_interval, double wavelength, GeoPos3D earth_point, int n_azimuth,
    const GeoPos3D* sensor_position, const Velocity3D* sensor_velocity)
{
    // binary search is used in finding the zero doppler time
    int lower_bound = 0;
    int upper_bound = n_azimuth - 1;
    auto lower_bound_freq = CalcDopplerFrequency(earth_point, sensor_position[lower_bound], sensor_velocity[lower_bound], wavelength);
    auto upper_bound_freq = CalcDopplerFrequency(earth_point, sensor_position[upper_bound], sensor_velocity[upper_bound], wavelength);

    if (std::abs(lower_bound_freq) < 1.0) {
        return lower_bound * line_time_interval;
    } else if (std::abs(upper_bound_freq) < 1.0) {
        return upper_bound * line_time_interval;
    } else if (lower_bound_freq * upper_bound_freq > 0.0) {
        return INVALID_AZ_TIME;
    }

    if constexpr (ZERO_DOPPLER_INTERCEPT) {
        double k = (upper_bound_freq - lower_bound_freq) / ((n_azimuth)*line_time_interval);
        double C = upper_bound_freq - k * n_azimuth * line_time_interval;
        double az_idx = -C / k;
        int az_idx_i = az_idx / line_time_interval;
        upper_bound = az_idx_i + 5;
        lower_bound = az_idx_i - 5;
    }

    // start binary search
    double mid_freq;
    while (upper_bound - lower_bound > 1) {
        const auto mid = (int)((lower_bound + upper_bound) / 2.0);
        mid_freq = sensor_velocity[mid].x * (earth_point.x - sensor_position[mid].x) + sensor_velocity[mid].y * (earth_point.y - sensor_position[mid].y) + sensor_velocity[mid].z * (earth_point.z - sensor_position[mid].z);

        if (mid_freq * lower_bound_freq > 0.0) {
            lower_bound = mid;
            lower_bound_freq = mid_freq;
        } else if (mid_freq * upper_bound_freq > 0.0) {
            upper_bound = mid;
            upper_bound_freq = mid_freq;
        } else if (mid_freq == 0.0) {
            return mid * line_time_interval;
        }
    }

    const auto y0 = lower_bound - lower_bound_freq * (upper_bound - lower_bound) / (upper_bound_freq - lower_bound_freq);

    // fmt::print(" bin search = {} , slope = {}, delta {}\n", y0 * line_time_interval, az_calc_idx, az_calc_idx - y0 * line_time_interval);
    return y0 * line_time_interval;
}

template <class T>
__device__ constexpr const T& min(const T& a, const T& b)
{
    return (b < a) ? b : a;
}

template <class T>
__device__ constexpr const T& max(const T& a, const T& b)
{
    return (a < b) ? b : a;
}

template <class T>
__device__ constexpr T blerp(
    T q11, T q21,
    T q12, T q22,
    T x, T y)
{
    // x and y are assumed to be in [0, 1]
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>);
    T bot = q11 + (q21 - q11) * x;
    T top = q12 + (q22 - q12) * x;

    return bot + (top - bot) * y;
}

[[maybe_unused]]
__device__ GeoPos3D GetPosition(double time, const OSV* vectors, int n_osv)
{
    int i0 = 0;
    int iN = n_osv - 1;
    if (iN > 8) {
        double t_first = vectors[0].tp;
        i0 = max<int>(0, (time - t_first) - 4);
        iN = min(iN, i0 + 8);
        // fmt::print("t = {} i0 iN = {} {}\n", time, i0, iN);
    }

    GeoPos3D result { 0, 0, 0 };
    for (int i = i0; i <= iN; ++i) {
        auto const orbI = vectors[i];

        double weight = 1;
        for (int j = i0; j <= iN; ++j) {
            if (j != i) {
                double const time2 = vectors[j].tp;

                weight *= (time - time2) / (orbI.tp - time2);
            }
        }
        result.x += weight * orbI.xp;
        result.y += weight * orbI.yp;
        result.z += weight * orbI.zp;
    }
    return result;
}

__global__ void RunTC(const float* RESTRICT data_in, const float* RESTRICT dem, Dimensions io_dims, Args args, float* RESTRICT data_out)
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

    const double az_time = GetEarthPointZeroDopplerTime(args.line_time_interval, args.wavelength, earth_point,
        io_dims.in_y_size, args.sat_pos_az, args.sat_vel_az);

    // const int out_idx = x + y * out_x_size;
    if (az_time == INVALID_AZ_TIME) {
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
    GeoPos3D sat_pos = { };
    sat_pos.x = pos0.x + interp_dt * (pos1.x - pos0.x);
    sat_pos.y = pos0.y + interp_dt * (pos1.y - pos0.y);
    sat_pos.z = pos0.z + interp_dt * (pos1.z - pos0.z);
#endif
    // range idx
    double dx = earth_point.x - sat_pos.x;
    double dy = earth_point.y - sat_pos.y;
    double dz = earth_point.z - sat_pos.z;
    double slant_range = sqrt(dx * dx + dy * dy + dz * dz);
    double rg_idx = (slant_range - args.slant_range_first_sample) * rg_m_to_pix;
    if (rg_idx < 0 || rg_idx >= in_x_size) {
        data_out[out_idx] = 0;
        return;
    }
    /*
    double slr_time = (2 * slant_range) / 299792458.0;
    if((y == (y_start + y_end) / 2) && (x == out_x_size / 2))

        fmt::print("az time = {} rg idx = {} slr = {}\n", az_time, slr_time, slant_range);
    */

#if INTERP_MODE == 0
    int rg_idx_nn = std::round(rg_idx);
    int az_idx_nn = std::round(az_idx);

    size_t in_idx = rg_idx_nn + in_x_size * az_idx_nn;
    data_out[out_idx] = data_in[in_idx];
#elif INTERP_MODE == 1
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

    Args args = { };
    args.lat_start = dem.gt[3];
    args.lon_start = dem.gt[0];
    args.pixel_spacing_x = dem.gt[1];
    args.pixel_spacing_y = dem.gt[5];

    args.no_data_value = dem.no_data_value;

    std::vector<GeoPos3D> pos;
    std::vector<Velocity3D> vel;
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

    DeviceBuffer<GeoPos3D> d_pos;
    d_pos.init(pos.size());
    d_pos.cpy_from_host(pos.data());

    DeviceBuffer<Velocity3D> d_vel;
    d_vel.init(vel.size());
    d_vel.cpy_from_host(vel.data());

    DeviceBuffer<OSV> d_osv;

    d_osv.init(osv.size());
    d_osv.cpy_from_host(osv.data());

    args.line_time_interval = sar_meta.line_time_interval;
    args.wavelength = sar_meta.wavelength;
    args.range_spacing = sar_meta.range_spacing;
    args.sat_pos_az = d_pos.data();
    args.sat_vel_az = d_vel.data();
    args.osv = d_osv.data();
    args.n_osv = d_osv.size();
    args.slant_range_first_sample = sar_meta.slant_range_first_sample;

    Dimensions io_dims = { };
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
