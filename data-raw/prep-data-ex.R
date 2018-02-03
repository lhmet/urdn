
area_file <- "/home/pqgfapergs1/Dropbox/ufsm/orientacoes/pos-doc/nelson/work/arpfunc/DirFluxo4_test/AreaAlta.rst"
dir_file <- "/home/pqgfapergs1/Dropbox/ufsm/orientacoes/pos-doc/nelson/work/arpfunc/DirFluxo4_test/DirAlta.rst"
area_hres <- raster::raster(area_file)
area_hres
dir_hres <- raster::raster(dir_file)
dir_hres
dir_hres[is.na(dir_hres)] <- -9999
area_hres[is.na(area_hres)] <- -9999
use_data(area_hres, dir_hres, overwrite = TRUE)

