.First.lib <- function(lib, pkg) {
  # Load dll:
  library.dynam("urdn", pkg, lib)
}
