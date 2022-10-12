.onLoad = function(libname, pkgname){
	# setting some options


	options("fixest_dict" = c())
	options("fixest_notes" = TRUE)
	options("fixest_print" = list(type = "table"))
	options("fixest_fl_authorized" = FALSE)

	setFixest_coefplot("all", reset = TRUE)
	setFixest_ssc()
	setFixest_etable()

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

.onAttach = function(libname, pkgname) {
    rlang::inform("FIXEST")
    rlang::inform(c("i" = "This is a fork of the original project created by Laurent Berge."))
    rlang::inform(c("i" = "Visit https://buymeacoffee.com/pacha if you decide to donate and contribute to improve this project."))
}


