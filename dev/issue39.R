devtools::load_all()

# --- re-building ‘exporting_tables.Rmd’ using rmarkdown

est_slopes = feols(Ozone ~ Solar.R + Wind | Day + Month[Temp], airquality)

# --- re-building ‘fixest_walkthrough.Rmd’ using rmarkdown
gravity_pois = fepois(Euros ~ log(dist_km) | Origin + Destination + Product + Year, trade)
print(gravity_pois)

# --- re-building ‘standard_errors.Rmd’ using rmarkdown
data(trade)
# OLS estimation
gravity = feols(log(Euros) ~ log(dist_km) | Destination + Origin + Product + Year, trade)
# Two-way clustered SEs
summary(gravity, vcov = "twoway")
# Two-way clustered SEs, without small sample correction
summary(gravity, vcov = "twoway", ssc = ssc(adj = FALSE, cluster.adj = FALSE))
