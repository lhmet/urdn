
area_file <- here::here("../../../../../Dropbox/documents/pauliceia/AreaAlta.rst")
dir_file <- here::here("../../../../../Dropbox/documents/pauliceia/DirAlta.rst")
area_hres <- raster::raster(area_file)
area_hres
dir_hres <- raster::raster(dir_file)
dir_hres
dir_hres[is.na(dir_hres)] <- -9999
area_hres[is.na(area_hres)] <- -9999
usethis::use_data(area_hres, dir_hres, overwrite = TRUE)

