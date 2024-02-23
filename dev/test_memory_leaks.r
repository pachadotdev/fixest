if (!require("devtools")) install.packages("devtools")
devtools::load_all()

# is_DT <- requireNamespace("data.table", quietly = TRUE)
# if (is_DT) library(data.table)

# logit, trade, quakes, did, smallsample, othervcovs, ivs, interactions, formulas, dotsquare
tests <- c(T, rep(F, 9))

# LOGIT ----

debug(cpp_demean)
fit <- feols(mpg ~ wt | cyl, mtcars)
fml <- mpg ~ wt | cyl
data <- mtcars
undebug(cpp_demean)

if (tests[1]) {
  fit <- fixest::feglm(am ~ wt | cyl, mtcars, family = binomial)
  s1 <- summary(fit)
}

# TRADE ----

if (tests[2]) {
  # if (!require("dplyr")) install.packages("dplyr")
  # library(dplyr)

  # trade_short <- Trade %>%
  #   filter(Year %in% c(2015, 2016)) %>%
  #   group_by(Year, Destination, Origin) %>%
  #   summarise(
  #     Euros = sum(Euros, na.rm = TRUE),
  #     dist_km = mean(dist_km, na.rm = TRUE)
  #   ) %>%
  #   ungroup() %>%
  #   as.data.frame()

  # use_data(trade_short, overwrite = TRUE)

  gravity_pois <- fepois(Euros ~ log(dist_km) | Origin + Destination + Year, trade_short)
  s1 <- summary(gravity_pois, vcov = "twoway")
  s2 <- summary(gravity_pois, vcov = ~Year)
  s3 <- summary(gravity_pois, cluster = "Year")
  s4 <- summary(gravity_pois, cluster = ~Year)
  s5 <- summary(gravity_pois, cluster = ~Year)

  gravity_simple <- fepois(Euros ~ log(dist_km), trade_short)
  s6 <- summary(gravity_simple, ~ Origin + Destination)

  gravity_pois_2 <- fepois(Euros ~ log(dist_km), trade_short, vcov = ~Year)

  gravity_ols <- feols(log(Euros) ~ log(dist_km) | Origin + Destination + Year, trade_short)

  gravity_negbin <- fenegbin(Euros ~ log(dist_km) | Origin + Destination + Year, trade_short)

  res_multi <- fepois(Euros ~ log(dist_km) | csw0(Year, Destination, Origin), trade_short)

  gravity_pois_fes <- fixef(gravity_pois)
  s7 <- summary(gravity_pois_fes)
}

if (tests[3]) {
  fit1 <- feols(depth ~ mag, quakes, "conley")
  fit2 <- feols(depth ~ mag, quakes, conley(200, distance = "spherical"))
  fit3 <- feols(depth ~ mag, quakes, vcov_conley(
    lat = "lat", lon = "long",
    cutoff = 200, distance = "spherical"
  ))
}

if (tests[4]) {
  est1 <- feols(y ~ x1, base_did)
  s1 <- summary(est1, newey ~ id + period)

  est2 <- feols(y ~ x1, base_did, panel.id = ~ id + period)
  s2 <- summary(est2, "newey_west")

  setFixest_estimation(panel.id = ~ id + period)
  est3 <- feols(y ~ x1, base_did)
  s3 <- summary(est3, "newey_west")
  s4 <- summary(est3, "cluster")

  est4 <- feols(y ~ x1 | period, base_did, "cluster")

  setFixest_estimation(reset = TRUE)

  est5 <- feols(y ~ x1 | period, base_did, "cluster")
  est6 <- feols(y ~ x1 | period, base_did, ~ id + period)
  est7 <- feols(y ~ x1 | period, base_did, NW(2) ~ id + period)
  est8 <- feols(y ~ x1 | period, base_did, vcov_NW("id", "period", lag = 2))
}

# SMALL SAMPLE CORRECTION ----

if (tests[5]) {
  est <- feols(y ~ x1 | id, base_did)
  est_up <- feols(y ~ x1 | id, base_did, ssc = ssc(fixef.K = "full"))
  est_down <- feols(y ~ x1 | id, base_did,
    ssc = ssc(adj = FALSE, cluster.adj = FALSE)
  )

  m1 <- feols(y ~ x1 | id, base_did, iid ~ ssc(adj = FALSE))
  m2 <- feols(y ~ x1 | id, base_did, hetero ~ ssc(adj = FALSE))
}

# OTHER VCOVS ----

if (tests[6]) {
  est <- feols(y ~ x1 | id, base_did)
  summ <- summary(est, vcov = sandwich::vcovHC, type = "HC1")
  fit <- feols(y ~ x1 | id, base_did, vcov = function(x) sandwich::vcovHC(x, type = "HC1"))
}

# IVS ----

if (tests[7]) {
  base <- iris
  names(base) <- c("y", "x1", "x_endo_1", "x_inst_1", "fe")
  set.seed(2)
  base$x_inst_2 <- 0.2 * base$y + 0.2 * base$x_endo_1 + rnorm(150, sd = 0.5)
  base$x_endo_2 <- 0.2 * base$y - 0.2 * base$x_inst_1 + rnorm(150, sd = 0.5)

  est_iv <- feols(y ~ x1 | x_endo_1 + x_endo_2 ~ x_inst_1 + x_inst_2, base)

  fs <- fitstat(est_iv, ~ ivf1 + ivwald1 + ivf2 + ivwald2, cluster = "fe")

  setFixest_print(fitstat = ~ . + ivwald2)
  est_iv

  est_iv_fe <- feols(
    y ~ x1 | fe | x_endo_1 + x_endo_2 ~ x_inst_1 + x_inst_2,
    base
  )

  s1 <- summary(est_iv_fe, stage = 1)
}

# INTERACTION TERMS ----

if (tests[8]) {
  base <- iris
  names(base) <- c("y", paste0("x", 1:3), "fe1")
  base$fe2 <- rep(letters[1:5], 30)
  est_comb <- feols(y ~ x1 | fe1^fe2, base)
  fes <- fixef(est_comb)[[1]]
  est_vs <- feols(y ~ x1 | fe1[x2], base)
  s_vs <- summary(fixef(est_vs))

  fit <- feols(Ozone ~ Solar.R + i(Month), airquality)

  res_i1 <- feols(Ozone ~ Solar.R + i(Month), airquality)
  res_i2 <- feols(Ozone ~ Solar.R + i(Month, ref = 8), airquality)
  res_i3 <- feols(Ozone ~ Solar.R + i(Month, keep = 5:6), airquality)

  est_did <- feols(y ~ x1 + i(period, treat, 5) | id + period, base_did)

  res_twfe <- feols(y ~ x1 + i(time_to_treatment, ref = c(-1, -1000)) |
    id + year, base_stagg)

  res_sa20 <- feols(y ~ x1 + sunab(year_treated, year) | id + year, base_stagg)

  s_sa20 <- summary(res_sa20, agg = "att")
}

# FORMULAS ----

if (tests[9]) {
  base <- iris
  names(base) <- c("y", "x1", "x2", "x3", "species")

  setFixest_fml(..ctrl = ~ poly(x2, 2) + poly(x3, 2))

  fit1 <- xpd(y ~ x1 + ..ctrl)

  vars <- c("x2", "x2^2", "x3")
  for (i in 1:3) {
    fit2 <- xpd(y ~ x1 + ..ctrl, ..ctrl = vars[1:i])
  }

  fit3 <- feols(y ~ x1 + ..ctrl, base)

  fit4 <- xpd(Armed.Forces ~ Population + regex("GNP|ployed"), data = longley)
}

# DOT SQUARE ----

if (tests[10]) {
  base <- setNames(iris, c("y", "x1", "x2", "x3", "species"))
  i <- 2:3
  z <- "i(species)"
  fit <- feols(y ~ x.[i] + .[z], base)

  i <- 1:3
  fit2 <- xpd(y ~ .["x.[i]_sq"])
}
