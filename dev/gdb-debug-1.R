data(airquality)
fixest2::feols(Ozone ~ Solar.R | csw(Month, Day), data = airquality)
