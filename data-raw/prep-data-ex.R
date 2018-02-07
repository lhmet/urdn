
area_file <- here::here("../../../../../Dropbox/documents/pauliceia/AreaAlta.rst")
dir_file <- here::here("../../../../../Dropbox/documents/pauliceia/DirAlta.rst")
area_hres <- raster::raster(area_file)
area_hres
dir_hres <- raster::raster(dir_file)
dir_hres
dir_hres[dir_hres < 0] <- NA
area_hres[area_hres < 0] <- NA
usethis::use_data(area_hres, dir_hres, overwrite = TRUE)

