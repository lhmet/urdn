#' Upscaling river drainage network
#'
#' Upscaling flow direction from fine- to coarse‐resolution river networks
#'
#' @param dirh high resolution flow direction grid.
#' @param areah high resolution drainage area grid.
#' @param factor integer. Aggregation factor expressed as number of cells in
#' each direction (horizontally and vertically). For example, fact=2 will result
#' in a new Raster* object with 2*2=4 times fewer cells.
#' @param area.thres area threshold (set equal to the cell area).
#' @param mufp minimum upstream flow path (set as one fifth of the cell size).
#' @param mascfile mask file
#'
#' @return list object with 2 elements:
##'  \describe{
##'   \item{hires}{RasterLayer with location of cells outlet}
##'   \item{lowres}{RasterStack with flow direction, outlet row, and outlet col}
##'   }
#' @export
#' @examples
#'area_hres
#'dir_hres
#'flowd_coarse_grid <- cotat_plus(dirh = dir_hres,
#'                                areah = area_hres,
#'                                factor = 10,
#'                                area.thres = 1.0,
#'                                mufp = 0.02)
#'# outlet locations (high resolution)
#'flowd_coarse_grid[[1]]
#'# flow direction, outlet row, outlet col
#'flowd_coarse_grid[[2]]
#'plot(flowd_coarse_grid[[2]])
#'flow_dir <- raster(flowd_coarse_grid[[2]], 1)
#'flow_dir

cotat_plus <- function(dirh,
                       areah,
                       factor,
                       area.thres,
                       mufp,
                       mascfile = "") {

  # # TO TEST
  # # if lib does not exist, create it
  # if (!file.exists("src/fluxdir.so") |
  #     !file.exists("libs/fluxdir.so")) {
  #   libf <- stringr::str_replace("src/fluxdir.so", "\\.so", "\\.f90")
  #   stopifnot(file.exists(libf))
  #   system(paste0("R CMD SHLIB ", libf))
  #   rm(libf)
  #   lib_path <- system.file(package = "urdn", "src/fluxdir.so")
  #   dyn.load(lib_path)
  #   stopifnot(is.loaded("fluxdir"))
  # }
  # #on.exit(dyn.unload(lib))

  nrowh <- nrow(dirh)
  ncolh <- ncol(dirh)

  if (!(nrowh %% factor == 0 & ncolh %% factor == 0)) {
    stop("high-resolution is not divisible by factor.")
  }

  nptsh <- nrowh * ncolh
  resh <- raster::xres(dirh)

  # Set low-resolution dimensions and parameters
  nrowl <- nrowh %/% factor
  ncoll <- ncolh %/% factor
  nptsl <- nrowl * ncoll
  resl <- factor * resh
  reshd <- resh / 1e5

  if (missing(area.thres)) {
    # default value
    area.thres <- resl / 1000.0
  }
  if (missing(mufp)) {
    # default value
    mufp <- (resh / 1000.0) / 5.0
  }
  pathlim <- floor(mufp / 100.0 / resh)

  # Load/create mask file
  # Load/create masc file
  if (missing(mascfile)) {
    masc <- raster::aggregate(dirh, factor)
    masc <- raster::setValues(masc, rep.int(0L, nptsl))
  } else {
    masc <- raster::raster(mascfile)
    if (raster::nrow(masc) != nrowl | raster::ncol(masc) != ncoll) {
      stop("masc dimensions don't match low-res data raster.")
    }
  }

  # Call the Fortran subroutine
  ans <- .Fortran(
    "fluxdir",
    as.integer(nrowh),
    as.integer(ncolh),
    as.integer(raster::getValues(dirh)),
    as.single(raster::getValues(areah)),
    outletm = integer(nptsh),
    as.integer(nrowl),
    as.integer(ncoll),
    as.integer(raster::getValues(masc)),
    dirlm = integer(nptsl),
    outrowm = integer(nptsl),
    outcolm = integer(nptsl),
    as.integer(factor),
    as.integer(factor),
    as.single(area.thres),
    as.integer(pathlim),
    PACKAGE = "urdn"
  )

  # Set output objects
  outlet <- raster::raster(dirh)
  dirl <- raster::raster(masc)
  outrow <- raster::raster(masc)
  outcol <- raster::raster(masc)

  outlet <- raster::setValues(outlet, ans$outletm)
  dirl <- raster::setValues(dirl, ans$dirlm)
  outrow <- raster::setValues(outrow, ans$outrowm)
  outcol <- raster::setValues(outrow, ans$outcolm)
  rm(ans)

  s <- raster::stack(dirl, outrow, outcol)
  names(s) <- c("flowdir", "outlet_row", "outlet_col")
  out_list <- c(hires = outlet, lowres = s)
  # a list with outlet
  return(out_list)
}