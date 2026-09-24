cmake_minimum_required(VERSION 3.25)
project(sartcpp)


include(FetchContent)

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

set(CMAKE_CXX_STANDARD 23)

add_executable(sartcpp main.cpp
        s1/s1_parser.cpp
        s1/s1_parser.hpp
        sar/sar_metadata.hpp
        s1/s1_raster_load.cpp
        s1/s1_raster_load.hpp
        util/memory_raster.hpp
        sar/dem.hpp
        sar/sar_geocode.cpp
        sar/sar_geocode.hpp
        util/math_utils.hpp
        util/proj_utils.hpp
        s1/s1_metadata.hpp
        s1/s1_deburst.cpp
        s1/s1_deburst.hpp
        util/gdal_util.hpp
        s1/s1_calibrate.cpp
        s1/s1_calibrate.hpp)

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


target_compile_definitions(sartcpp PRIVATE BOOST_DATE_TIME_NO_LIB)

find_package(Boost CONFIG REQUIRED COMPONENTS date_time)

target_link_libraries(sartcpp PRIVATE fmt::fmt pugixml gdal pthread Boost::date_time)
