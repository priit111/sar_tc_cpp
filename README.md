# SAR TC CPP

Synthetic aperture radar terrain correction implementation in C++(and Cuda)


Implements the following Sentinel1 processing scenario using SNAP terminology:

S1 Product read -> LVL1 SLC Split -> Calibrate(Beta) -> Deburst -> Range Doppler Terrain Correction -> Output GTIFF


* C++(CPU) implementation with no legacy mistakes in  with an eye for performance baseline comparison. 
* Implements the same processing pipeline with Cuda for a good CPU vs GPU performance comparison
* Minimal barebones implementation with easy build setup, the S1 metadata parser is minimal. Error checking is minimal as well. Does not intend to cover a full-fledged processor's capabilities.
* Only supports IW mode 1 swath 1 polarization only for ease of implementation
* Intended as a proof of concept processor for discussion about SAR processor design and performance with a concrete implementation instead of a common pipeline instead of abstract ideas
* Example discussion points when discussing with other developers:
  * C vs C++ vs rust vs Java vs python
  * Memory usage
  * Run time performance
  * GPU support - if, when and how to implement 



## Usage

./sartcpp {PATH_TO_S1_FOLDER_ROOT} {PATH_TO_DEM} {polarization} {swath} {output_path_name}

sartcpp ./S1A_IW_SLC__1SDV_20260522T155647_20260522T155714_064632_082444_F771.SAFE ./eesti_soome.tif vv iw1 /tmp/tc.tif


cuda version is the same with the binary name being sartcuda


## DEM usage

Idea copied from sarsen(https://github.com/bopen/sarsen), the processor simply uses the input DEM for output file generation. Simplifies the processor implementation(no DEM interpolation).

### DEM generation
Todo insert example script



## Building

```
mkdir build
cd build
cmake ..
make
```

Requires the following libraries:
```
GDAL
boost - plan is to remove it in the future
cuda - if building the GPU variant)
```

Everything else is done via fetchcontent.


# TODO / Ideas for the future
## RTC 
Implement either Small, Shiroma or both. SNAP calls in Terrain Flattening Op. More complicated and has again performance and memory usage implications.

## Interferometry
Again know to be slow/memory hungry, could be interesting to experiment in the future.

## Further optimizations

### algorithmic optimizations
TC implementation is same as SNAP(reference Guide-to-Sentinel-1-Geocoding.pdf) with a few tricks, but the overall main change is the approach to tiling and overall memory usage. The details here can be optimized and the memory access patterns and how the threading is done seem to have a major impact on TC run time. Additionally, an idea is to reimplement the TC zero Doppler condition finding with via alternate methods. Hard to tell if it faster or not.

### CPU
Main question seems to be the TC implementation threading and memory access patterns

### GPU optimizations
Datacenter vs consumer(with poor double FLOPS). Probably the improvement thing here is to experiment with nvTiff and improving the allocations and H->D and D->H transfers.

