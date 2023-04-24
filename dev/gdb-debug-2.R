data(airquality)
fixest2::feols(Ozone ~ Solar.R + sw0(Wind + Temp) | csw(Month, Day), airquality, cluster = ~Day)
