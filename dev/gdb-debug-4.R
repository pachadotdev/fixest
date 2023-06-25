library(fixest2)

load_all()

gravity_pois = fepois(Euros ~ log(dist_km) | Origin + Destination + Product + Year, trade)
fixedEffects = fixef(gravity_pois)

# fixef.fixest function L882
object = gravity_pois
notes = getFixest_notes()
sorted = TRUE
nthreads = getFixest_nthreads()
fixef.tol = 1e-5
fixef.iter = 10000