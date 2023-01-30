.onLoad <- function(libname, pkgname) {
  # setting some options


  options("fixest_dict" = c())
  options("fixest_notes" = TRUE)
  options("fixest_print" = list(type = "table"))
  options("fixest_fl_authorized" = FALSE)

  # setFixest_coefplot("all", reset = TRUE) # TODO: remove coefplot
  setFixest_ssc()
  # setFixest_etable() # TODO: remove etable

  # # To include later
  # cpp_setup_fork_presence()

  # nthreads
  setFixest_nthreads()

  # Setup of builtin VCOVs
  vcov_setup()

  # Aliases must come after the VCOV
  create_aliases()

  # To circumvent a peculiar behavior from pkgdown
  fix_pkgwdown_path()

  invisible()
}

.onAttach <- function(libname, pkgname) {
    packageStartupMessage("=== fixest2 ===\nTHIS IS A WORK IN PROGRESS!\nDON'T USE THIS EXCEPT FOR TESTING\nThis is a fork of the original fixest project created by Laurent Berge.\nVisit https://buymeacoffee.com/pacha if you decide to donate and contribute to improve this project.")
}
