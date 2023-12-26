# library(fixest2)
devtools::load_all()

# REPRODUCIBLE ERROR ---

# gravity_pois = fepois(Euros ~ log(dist_km) | Origin + Destination + Product + Year, trade)
# fixedEffects = fixef(gravity_pois)
#  *** caught segfault ***
# address 0x55cf901cb924, cause 'memory not mapped'

# Traceback:
#  1: .Call(`_fixest2_cpp_get_fe_gnl_`, as.integer(Q), as.integer(N),     sumFE, dumMat, as.integer(cluster_sizes), obsCluster)
#  2: cpp_get_fe_gnl(Q, N, S, dumMat, nbCluster, orderCluster)
#  3: fixef.fixest(gravity_pois)
#  4: fixef(gravity_pois)

# Possible actions:
# 1: abort (with core dump, if enabled)
# 2: normal R exit
# 3: exit R without saving workspace
# 4: exit R saving workspace

# NOW WE SEE IT LINE BY LINE ----

# the problem is in fixef.fixest function L882 Methods.R
# https://github.com/pachadotdev/fixest2/blob/cpp11_wip/R/Methods.R#L882
# that function calls https://github.com/pachadotdev/fixest2/blob/cpp11_wip/src/05_01_misc_helpers.cpp#L540

# this is the same as to run fixef.fixest line by line
# lines 924-1079 that do not apply for this case

gravity_pois <- fepois(Euros ~ log(dist_km) | Origin + Destination + Product + Year, trade)

object <- gravity_pois
notes <- getFixest_notes()
sorted <- TRUE
nthreads <- getFixest_nthreads()
fixef.tol <- 1e-5
fixef.iter <- 10000
S <- object$sumFE
family <- object$family
fixef_names <- object$fixef_vars
fixef_id <- object$fixef_id
Q <- length(fixef_id)
N <- length(S)
id_dummies_vect <- list()
for (i in 1:Q) id_dummies_vect[[i]] <- as.vector(fixef_id[[i]])
is_ref_approx <- FALSE
isSlope <- FALSE
dumMat <- matrix(unlist(id_dummies_vect), N, Q) - 1
orderCluster <- matrix(unlist(lapply(id_dummies_vect, order)), N, Q) - 1
nbCluster <- sapply(fixef_id, max)

input <- list(Q = Q, N = N, S = S, dumMat = dumMat, nbCluster = nbCluster, orderCluster = orderCluster)
saveRDS(input, "dev/input.rds")
