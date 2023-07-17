# library(fixest2)
devtools::load_all()

# REPRODUCIBLE ERROR ---

# gravity_pois = fepois(Euros ~ log(dist_km) | Origin + Destination + Product + Year, trade)
# fixedEffects = fixef(gravity_pois)
# Error: Invalid input type, expected 'integer' actual 'double'

# NOW WE SEE IT LINE BY LINE ----

# the problem is in fixef.fixest function L882 Methods.R

# this is the same as to run fixef.fixest line by lines 924-1079 that do not
# apply for this case

gravity_pois = fepois(Euros ~ log(dist_km) | Origin + Destination + Product + Year, trade)

object <- gravity_pois
notes <- getFixest_notes()
sorted <- TRUE
nthreads <- getFixest_nthreads()
fixef.tol <- 1e-5
fixef.iter <- 10000

check_arg(notes, sorted, "logical scalar")

check_value(fixef.tol, "numeric scalar GT{0} LT{1}")
check_value(fixef.iter, "strict integer scalar GT{0}")

# Preliminary stuff
S <- object$sumFE

family <- object$family
fixef_names <- object$fixef_vars

fixef_id <- object$fixef_id

Q <- length(fixef_id)
N <- length(S)

# either (we need to clean its attributes for unlist to be efficient)
id_dummies_vect <- list()
for (i in 1:Q) id_dummies_vect[[i]] <- as.vector(fixef_id[[i]])

is_ref_approx <- FALSE
isSlope <- FALSE

# We apply a cpp script to handle complicated cases (and we don't know beforehand if the input is one)

dumMat <- matrix(unlist(id_dummies_vect), N, Q) - 1
orderCluster <- matrix(unlist(lapply(id_dummies_vect, order)), N, Q) - 1

nbCluster <- sapply(fixef_id, max)

print(paste("Q", class(Q)))
print(paste("N", class(N)))
print(paste("S", class(S)))
print(paste("dumMat", class(dumMat)))
print(paste("nbCluster", class(nbCluster)))
print(paste("orderCluster", class(orderCluster)))

print(head(dumMat))
print(head(orderCluster))

# ERROR 1 ----

fixef_values <- cpp_get_fe_gnl(Q, N, S, dumMat, nbCluster, orderCluster)

# OK L570OK L5954
# 0
# OK L612OK L6321

#  *** caught segfault ***
# address 0x5623b3cd0df0, cause 'memory not mapped'

# Traceback:
#  1: .Call(`_fixest2_cpp_get_fe_gnl_`, as.integer(Q), as.integer(N),     sumFE, as.integer(dumMat), as.integer(cluster_sizes), as.integer(obsCluster))
#  2: cpp_get_fe_gnl(Q, N, S, dumMat, nbCluster, orderCluster)

# Possible actions:
# 1: abort (with core dump, if enabled)
# 2: normal R exit
# 3: exit R without saving workspace
# 4: exit R saving workspace

# ERROR 2 ----

storage.mode(dumMat) <- "integer"
storage.mode(orderCluster) <- "integer"
fixef_values <- cpp_get_fe_gnl(Q, N, S, dumMat, nbCluster, orderCluster)

# OK L570OK L5954
# 0
# OK L612OK L6321

#  *** caught segfault ***
# address 0x56144c29706c, cause 'memory not mapped'

# Traceback:
#  1: .Call(`_fixest2_cpp_get_fe_gnl_`, as.integer(Q), as.integer(N),     sumFE, as.integer(dumMat), as.integer(cluster_sizes), as.integer(obsCluster))
#  2: cpp_get_fe_gnl(Q, N, S, dumMat, nbCluster, orderCluster)

# Possible actions:
# 1: abort (with core dump, if enabled)
# 2: normal R exit
# 3: exit R without saving workspace
# 4: exit R saving workspace
