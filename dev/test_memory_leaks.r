if (!require("devtools")) install.packages("devtools")
devtools::load_all()

# Simple example using trade data ----

data(trade)

# subset to 2015-2016
trade <- subset(trade, Year %in% c(2015, 2016))

gravity_ols <- feols(log(Euros) ~ log(dist_km) | Origin + Destination +
  Product + Year, trade)

gravity_pois <- fepois(Euros ~ log(dist_km) | Origin + Destination + Product +
  Year, trade)

gravity_negbin <- fenegbin(Euros ~ log(dist_km) | Origin + Destination +
  Product + Year, trade)

summary(gravity_ols, vcov = "twoway")
summary(gravity_ols, vcov = ~Product)
summary(gravity_ols, cluster = "Product")
summary(gravity_ols, cluster = ~Product)

summary(gravity_pois, vcov = "twoway")
summary(gravity_pois, vcov = ~Product)
summary(gravity_pois, cluster = "Product")
summary(gravity_pois, cluster = ~Product)

summary(gravity_negbin, vcov = "twoway")
summary(gravity_negbin, vcov = ~Product)
summary(gravity_negbin, cluster = "Product")
summary(gravity_negbin, cluster = ~Product)

fixef(gravity_ols)
fixef(gravity_pois)
fixef(gravity_negbin)

# Simple example using iris + factor partition ----

nm <- names(iris)

fit1 <- feols(.[nm[1]] ~ .[nm[2:4]], iris, fsplit = ~Species)
fit2 <- fepois(.[nm[1]] ~ .[nm[2:4]], iris, fsplit = ~Species)
fit3 <- fenegbin(.[nm[1]] ~ .[nm[2:4]], iris, fsplit = ~Species)
