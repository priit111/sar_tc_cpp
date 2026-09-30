// SPDX-License-Identifier: GPL-3.0-or-later
// https://github.com/priit111/sar_tc_cpp

#include "s1_parser.hpp"

#include <charconv>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <gdal/gdal_priv.h>
#include <iostream>

#include "pugixml.hpp"

#include "../util/proj_utils.hpp"
#include "s1_metadata.hpp"

#include "../sar/sar_metadata.hpp"

// void load_annot(std::string)

using namespace s1;

namespace {

template <class T>
std::vector<T> spaced_str_arr_to(std::string_view src)
{
    std::vector<std::string_view> tokens;

    size_t start_idx = 0;
    while (1) {
        auto next = src.find(" ", start_idx);
        if (next == std::string_view::npos) {
            tokens.push_back(src.substr(start_idx, src.size() - start_idx));
            break;
        }
        tokens.push_back(src.substr(start_idx, next - start_idx));
        start_idx = next + 1;
    }

    std::vector<T> ret;
    for (std::string_view e : tokens) {
        /*
         libstdc++ float support...
        T val;
        auto res = std::from_chars(e.begin(), e.end(), val);
        if (res.ec != std::errc { }) {
            SARTCPP_ASSERT(false);
        }
        */
        T val = std::strtod(e.begin(), nullptr);
        ret.push_back(val);
    }
    return ret;
}

double get_double(const pugi::xml_node& n)
{
    return n.text().as_double();
}

int get_int(const pugi::xml_node& n)
{
    return n.text().as_int();
}

std::string get_str(const pugi::xml_node& n)
{
    return n.text().as_string();
}

void load_cal(std::string path, S1Metadata& s1_meta)
{

    pugi::xml_document doc;
    auto res = doc.load_file(path.c_str());

    // doc.print(std::cout);

    auto cvl = doc.child("calibration").child("calibrationVectorList");

    for (const auto& el : cvl.children("calibrationVector")) {
        CalibMeta cm = { };
        std::string az_str = get_str(el.child("azimuthTime"));

        cm.az_time = parse_abs_time(az_str);
        cm.beta = spaced_str_arr_to<float>(get_str(el.child("betaNought")));
        cm.sigma = spaced_str_arr_to<float>(get_str(el.child("sigmaNought")));
        cm.gamma = spaced_str_arr_to<float>(get_str(el.child("gamma")));
        cm.pixel = spaced_str_arr_to<int>(get_str(el.child("pixel")));

        SARTCPP_ASSERT(cm.beta.size() == cm.gamma.size());
        SARTCPP_ASSERT(cm.beta.size() == cm.sigma.size());
        SARTCPP_ASSERT(cm.beta.size() == cm.pixel.size());
        // el.print(std::cout);
        s1_meta.calib.push_back(std::move(cm));
    }

    // cvl.print(std::cout);

    if (!res) {
        fmt::print("ERROR!\n");
    }
}

void load_annot(std::string path, SARMetadata& sar_meta, S1Metadata& s1_meta)
{
    pugi::xml_document doc;
    auto res = doc.load_file(path.c_str());
    if (!res) {
        fmt::print("ERROR parsing = {}!\n", path);
        SARTCPP_ASSERT(false);
    }
    // doc.print(std::cout);
    // doc.child("product").child("generalAnnotation").print(std::cout);

    auto ol = doc.child("product").child("generalAnnotation").child("orbitList");
    auto ii = doc.child("product").child("imageAnnotation").child("imageInformation");

    sar_meta.range_size = get_int(ii.child("numberOfSamples"));
    sar_meta.azimuth_size = get_int(ii.child("numberOfLines"));

    double slrt = get_double(ii.child("slantRangeTime"));

    const double c = 299792458.0;
    sar_meta.slant_range_first_sample = slrt * 0.5 * c;

    sar_meta.range_spacing = get_double(ii.child("rangePixelSpacing"));
    sar_meta.azimuth_spacing = get_double(ii.child("azimuthPixelSpacing"));
    sar_meta.line_time_interval = get_double(ii.child("azimuthTimeInterval"));
    double azimuth_spacing = get_double(ii.child("azimuthPixelSpacing"));
    std::string s = ii.child("productFirstLineUtcTime").text().as_string();

    auto first_line_time = parse_abs_time(s);

    sar_meta.frequency = get_double(doc.child("product").child("generalAnnotation").child("productInformation").child("radarFrequency"));
    sar_meta.wavelength = 299792458.0 / sar_meta.frequency;
    // fmt::print("wl = {}\n", sar_meta.wavelength);

    sar_meta.first_line_time = first_line_time;

    std::vector<OSV> osv_vec;
    for (auto osv_xml : ol.children()) {
        auto tp = parse_abs_time(osv_xml.child("time").text().as_string());

        OSV osv = { };

        osv.tp = (tp.timestamp_us - first_line_time.timestamp_us) / 1e6;
        {
            auto pos = osv_xml.child("position");
            osv.xp = get_double(pos.child("x"));
            osv.yp = get_double(pos.child("y"));
            osv.zp = get_double(pos.child("z"));
        }

        {
            auto vel = osv_xml.child("velocity");
            osv.xv = get_double(vel.child("x"));
            osv.yv = get_double(vel.child("y"));
            osv.zv = get_double(vel.child("z"));
        }
        osv_vec.push_back(osv);

        //fmt::print("[{}],[{} {} {}] [{} {} {}]\n", osv.tp, osv.xp, osv.yp, osv.zp, osv.xv, osv.yv, osv.zv);
    }
    sar_meta.osv = osv_vec;
    {
        auto st = doc.child("product").child("swathTiming");

        int lines_per_burst = get_int(st.child("linesPerBurst"));
        s1_meta.lines_per_burst = lines_per_burst;
        // fmt::print("lines_per_burst = [{}]\n", lines_per_burst);

        auto bl = st.child("burstList");

        for (auto e : bl.children("burst")) {
            // fmt::print("{}\n", e.value());

            std::string az_time = get_str(e.child("azimuthTime"));
            std::string first_valid_list = get_str(e.child("firstValidSample"));

            std::vector<int> vals = spaced_str_arr_to<int>(first_valid_list);
            BurstMeta bm = { };
            bm.az_time = parse_abs_time(az_time);
            bm.first_valid_sample = std::move(vals);
            s1_meta.bursts.push_back(std::move(bm));

            // fmt::print("AZ TIME = {}\nlist = {}\n", az_time, first_valid_list);
            // fmt::print("list size = {}\n", s1_meta.bursts.back().first_valid_sample.size());
        }
    }

    std::string geobox_str;
    {
        double min_lat = 100e3;
        double max_lat = -100e3;
        double min_lon = 100e3;
        double max_lon = -100e3;
        auto glgpl = doc.child("product").child("geolocationGrid").child("geolocationGridPointList");
        for (const auto& el : glgpl.children("geolocationGridPoint")) {
            //el.print(std::cout);
            double lat = get_double(el.child("latitude"));
            double lon = get_double(el.child("longitude"));
            GeoLocationGridPoint glgp = {};
            glgp.az_time = get_str(el.child("azimuthTime"));
            glgp.slrt = get_double(el.child("slantRangeTime"));
            glgp.line = get_int(el.child("line"));
            glgp.pixel = get_int(el.child("pixel"));
            glgp.lat = get_double(el.child("latitude"));
            glgp.lon = get_double(el.child("longitude"));
            glgp.height = get_double(el.child("height"));
            glgp.incidence = get_double(el.child("incidenceAngle"));
            glgp.elevation = get_double(el.child("elevationAngle"));
            s1_meta.geogrid_points[glgp.line].push_back(glgp);
            min_lat = std::min(lat, min_lat);
            max_lat = std::max(lat, max_lat);
            min_lon = std::min(lon, min_lon);
            max_lon = std::max(lon, max_lon);
        }

        /*
        for (const auto& e :  s1_meta.geogrid_points) {
            fmt::print("line = {}\n", e.first);
            for (const auto& el : e.second) {
                fmt::print("pix = {} inc = {}\n", el.pixel, el.incidence);
            }
        } */
        auto it = s1_meta.geogrid_points.begin();
        sar_meta.incidence_angle_begin = it->second.front().incidence;
        sar_meta.incidence_angle_end = it->second.back().incidence;
        geobox_str = fmt::format("Geobox = ({} {}) , ({} {})\n", min_lat, max_lat, min_lon, max_lon);
    }
    fmt::print("\n~~~ Metadata summary begin ~~~\n\n");
    fmt::print("rg spacing = {}\naz spacing = {}\n", sar_meta.range_spacing, sar_meta.azimuth_spacing);
    fmt::print("first line time = {}\n", abstime_to_str(sar_meta.first_line_time));
    fmt::print("center frequency = {}\n", sar_meta.frequency);
    fmt::print("wavelength = {}\n", sar_meta.wavelength);
    fmt::print("raster size = ({} , {})\n", sar_meta.range_size, sar_meta.azimuth_size);
    fmt::print("{}", geobox_str);

    {

        SARTCPP_ASSERT_MSG(!osv_vec.empty(), "OSV parse fail");
        fmt::print("middle OSV:\n");
        OSV osv = osv_vec[osv_vec.size()/2];
        fmt::print("[{}],[{} {} {}] [{} {} {}]\n", osv.tp, osv.xp, osv.yp, osv.zp, osv.xv, osv.yv, osv.zv);

    }

    fmt::print("\n~~~ Metadata summary end ~~~\n\n");
}
}

namespace s1 {
bool parse(std::string s1_dir_path, std::string pol, std::string swath, SARMetadata& sar_meta, S1Metadata& s1_meta)
{
    TimeBlock tb("metadata load");
    for (const auto& dirEntry : std::filesystem::recursive_directory_iterator(s1_dir_path)) {
        if (dirEntry.is_regular_file()) {

            std::string p = dirEntry.path().string();
            // fmt::print("path = {}\n", p);
            if (p.find(pol) != std::string::npos && p.find(swath) != std::string::npos) {
                if (p.find("rfi") != std::string::npos) {
                    continue;
                }
                if (p.find("noise") != std::string::npos) {
                    continue;
                }
                if (p.find("calibration") != std::string::npos) {
                    load_cal(p, s1_meta);
                    continue;
                }
                if (p.find(".tif") != std::string::npos) {
                    sar_meta.raster_path = p;
                    continue;
                }
                if (p.find(".xml") != std::string::npos) {
                    load_annot(p, sar_meta, s1_meta);
                }
                // fmt::print("{}\n", dirEntry.path().string());
            }
        }
    }

    //for now just assert...
    SARTCPP_ASSERT_MSG(!sar_meta.osv.empty(), "Metadata parsing problems...");
    SARTCPP_ASSERT_MSG(!s1_meta.bursts.empty(), "Metadata parsing problems...");
    SARTCPP_ASSERT_MSG(!s1_meta.calib.empty(), "Metadata parsing problems...");


    return true;
}
}