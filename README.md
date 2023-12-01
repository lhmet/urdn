# urdn

{*`urdn`*} is a R package for Upscaling River Drainage Networks.

The main function `cotat_plus()` is a wrapper for the algorithm coded in Fortran, known as Cell Outlet Tracing with an Area Threshold (COTAT+) (Paz et al. 2006) for upscaling flow direction from fine- to coarse‐resolution grid.

```{r}
library(urdn)
# raster of drainage area in high resolution
area_hres
# raster of flow direction in high resolution
dir_hres
# function usage
upscaling_out <- cotat_plus(flowdir.fres = dir_hres,
                           draina.fres = area_hres,
                           factor = 10 # factor to increase the spatial resolution 
)
# outlet locations (high resolution)
upscaling_out[[1]]
#flow direction, outlet row, outlet col
upscaling_out[[2]]
plot(upscaling_out[[2]])
flow_dir <- raster(upscaling_out[[2]], 1)
flow_dir
```

