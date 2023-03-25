devtools::load_all()
fenegbin(Euros ~ log(dist_km) | Origin + Destination + Product + Year, trade)
