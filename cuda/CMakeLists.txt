cmake_minimum_required(VERSION 3.25)
project(sartcuda CUDA CXX)


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

add_executable(sartcuda main_gpu.cpp
        ../s1/s1_parser.cpp
        ../s1/s1_parser.hpp
        ../sar/sar_metadata.hpp
        ../s1/s1_raster_load.cpp
        ../s1/s1_raster_load.hpp
        ../util/memory_raster.hpp
        ../sar/dem.hpp
        ./sar/sar_geocode_gpu.cu
        ../sar/sar_geocode.hpp
        ../util/math_utils.hpp
        ../util/proj_utils.hpp
        ../s1/s1_metadata.hpp
        ./s1/s1_deburst_gpu.cu
        ../s1/s1_deburst.hpp
        ../util/gdal_util.hpp
        ./s1/s1_calibrate_gpu.cu
        ../s1/s1_calibrate.hpp
)

target_include_directories(sartcuda PRIVATE "../")


find_package(CUDAToolkit 11.2 REQUIRED)
INCLUDE_DIRECTORIES(${CUDAToolkit_INCLUDE_DIRS})


find_package(Boost CONFIG REQUIRED COMPONENTS date_time)

target_link_libraries(sartcuda PRIVATE fmt::fmt pugixml gdal pthread Boost::date_time CUDA::cudart)
