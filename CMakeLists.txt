cmake_minimum_required(VERSION 3.20)
project(sartcpp)



set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

set(FETCHCONTENT_QUIET OFF)

FetchContent_Declare(
        pugixml
        GIT_REPOSITORY https://github.com/zeux/pugixml
)

FetchContent_MakeAvailable(pugixml)

include(FetchContent)

FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt
        GIT_TAG e69e5f977d458f2650bb346dadf2ad30c5320281
)
FetchContent_MakeAvailable(fmt)

FetchContent_Declare(
        Eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE)
FetchContent_MakeAvailable(Eigen)


add_executable(sartcpp main.cpp
        s1/s1_parser.cpp
        s1/s1_parser.hpp
        sar/sar_metadata.hpp
        s1/s1_raster_load.cpp
        s1/s1_raster_load.hpp
        util/memory_raster.hpp
        sar/dem.hpp
        sar/terrain_correction.cpp
        sar/terrain_correction.hpp
        util/math_utils.hpp
        util/proj_utils.hpp
        s1/s1_metadata.hpp
        s1/s1_deburst.cpp
        s1/s1_deburst.hpp
        util/gdal_util.hpp
        util/math_utils.hpp
        util/math_utils.cpp
        s1/s1_calibrate.cpp
        s1/s1_calibrate.hpp
        sar/orbit.hpp
        sar/time_type.hpp
        sar/multilook.cpp
        sar/multilook.hpp
        sar/rtc.cpp
        sar/rtc.hpp
        sar/sar_geo.hpp)

if(NOT CMAKE_BUILD_TYPE)
    message("Defaulting to release build...")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type." FORCE)
endif()


if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native")
    message("\n\n\n~~~Optimized build with -Ofast & -march=native!!!~~~\n\n\n")
else ()
    message("\n\n\n~~~Debug build!!!~~~\n\n\n")
endif ()



target_link_libraries(sartcpp PRIVATE fmt::fmt pugixml gdal pthread Eigen3::Eigen)
