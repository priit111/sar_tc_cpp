# SAR TC CPP

Synthetic aperture radar terrain correction implementation in C++(and Cuda)


Implements the following Sentinel1 processing scenario using SNAP terminology:

S1 Product read -> LVL1 SLC Split -> Calibrate(Beta) -> Deburst -> Range Doppler Terrain Correction -> Output GTIFF


* C++(CPU) implementation with no legacy mistakes in with an eye for performance baseline. 
* Implements the same processing pipeline with Cuda for a good CPU vs GPU performance comparison
* Minimal barebones implementation with easy build setup, the S1 metadata parser is minimal. Error checking is minimal as well. Does not intend to cover a full-fledged processor's capabilities.
* Only supports IW mode SLC 1 swath 1 polarization only due to ease of implementation
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

Idea copied from sarsen(https://github.com/bopen/sarsen), the processor simply uses the input DEM for output file generation. Simplifies the processor implementation(no DEM interpolation).  At the moment assumes epsg 4326 without checks. Proper projection handling is TODO.

### DEM generation

Example DEM creation via [pyroSAR](https://pypi.org/project/pyroSAR/)

```python
import os
from pyroSAR.auxdata import dem_autoload, dem_create
from spatialist.vector import bbox
from spatialist.auxil import utm_autodetect

# writeable proj lib path
os.environ["PROJ_LIB"] = ("./proj/")

# geobox of an example S1 one swath
extent = {'xmin': 23.14, 'xmax': 25.37, 'ymin': 59.1, 'ymax': 60.77}

vrt = 'my_test.vrt'
dem = 'my_test.tif'

# ~10m spacing at 60 latitude
tr = (18e-5,9e-5)


# download all needed tiles and mosaic them in a VRT (GDAL virtual file format)
with bbox(coordinates=extent, crs=4326) as geom:
    utm_epsg = utm_autodetect(geom, 'epsg')  # auto-detect the UTM zone
    dem_autoload(geometries=[geom], vrt=vrt, demType='Copernicus 30m Global DEM')

# create the EPSG:4326 DEM
# let's directly use the Copernicus DEM. This one uses EGM2008 geoid as vertical datum
dem_create(src=vrt, dst=dem, geoid_convert=True, geoid='EGM2008', pbar=True,
           tr=tr, nodata=-32768)


```



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
cuda - if building the GPU variant
```

Everything else is done via fetchcontent.


# TODO / Ideas for the future
## RTC 
WIP CPU D.Small RTC implemented, cuda TODO, Shiroma AP optimization TODO.

## Further optimizations

### algorithmic optimizations
Main variation in TC seems to be the zero Doppler finding function. Implemented(proj_conf.hpp): Bisection(from [SNAP](https://github.com/senbox-org/microwave-toolbox/blob/177cec6f66846e62b3d9471c7386311f7ecdcbee/sar-commons/src/main/java/eu/esa/sar/commons/SARGeocoding.java#L170)) and Newton approximation. 

### CPU
Main question seems to be the TC implementation threading and memory access patterns, 

### GPU optimizations
Datacenter vs consumer(with poor double FLOPS). Probably the improvement thing here is to experiment with nvTiff and improving the allocations and H->D and D->H transfers.

