devtools::load_all()
gravity_negbin = fenegbin(Euros ~ log(dist_km) | Origin + Destination + Product + Year, trade)
