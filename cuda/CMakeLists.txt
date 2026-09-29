cmake_minimum_required(VERSION 3.25)
project(sartcuda CUDA CXX)


set(CMAKE_CXX_STANDARD 17)



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

FetchContent_Declare(
        Eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE)
FetchContent_MakeAvailable(Eigen)

add_executable(sartcuda main_gpu.cpp
        ../s1/s1_parser.cpp
        ../s1/s1_parser.hpp
        ../sar/sar_metadata.hpp
        ../s1/s1_raster_load.cpp
        ../s1/s1_raster_load.hpp
        ../util/memory_raster.hpp
        ../sar/dem.hpp
        sar/terrain_correction_gpu.cu
        ../sar/terrain_correction.hpp
        ../util/math_utils.hpp
        ../util/math_utils.cpp
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

#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-Wall")
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-Werror")
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-Wextra")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr")


target_link_libraries(sartcuda PRIVATE fmt::fmt pugixml gdal pthread Eigen3::Eigen CUDA::cudart)
