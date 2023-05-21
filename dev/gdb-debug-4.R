library(fixest2)
gravity_pois = fepois(Euros ~ log(dist_km) | Origin + Destination + Product + Year, trade)
fixedEffects = fixef(gravity_pois)
