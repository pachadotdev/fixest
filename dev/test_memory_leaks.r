if (!require("devtools")) install.packages("devtools")
devtools::load_all()

# is_DT <- requireNamespace("data.table", quietly = TRUE)
# if (is_DT) library(data.table)

eval_trade <- F
eval_quakes <- F
eval_did <- F
eval_smallsample <- F
eval_othervcovs <- F
eval_ivs <- F
eval_interactions <- F # error with feols(Ozone ~ Solar.R + i(Month), airquality)
eval_formulas <- T

# TRADE ----

if (eval_trade) {
  # subset Year to 2015/16
  trade <- trade[Year %in% c(2015, 2016)]

  gravity_pois <- fepois(Euros ~ log(dist_km) | Origin + Destination + Product + Year, trade)
  s1 <- summary(gravity_pois, vcov = "twoway")
  s2 <- summary(gravity_pois, vcov = ~Product)
  s3 <- summary(gravity_pois, cluster = "Product")
  s4 <- summary(gravity_pois, cluster = ~Product)
  s5 <- summary(gravity_pois, cluster = ~Product)

  gravity_simple <- fepois(Euros ~ log(dist_km), trade)
  s6 <- summary(gravity_simple, ~ Origin + Destination)

  gravity_pois_2 <- fepois(Euros ~ log(dist_km), trade, vcov = ~Product)

  gravity_ols <- feols(log(Euros) ~ log(dist_km) | Origin + Destination + Product + Year, trade)

  gravity_negbin <- fenegbin(Euros ~ log(dist_km) | Origin + Destination + Product + Year, trade)

  tab1 <- etable(gravity_pois, gravity_negbin, gravity_ols,
    vcov = "twoway", headers = c("Poisson", "Negative Binomial", "Gaussian")
  )

  tab2 <- etable(gravity_pois, gravity_ols, vcov = "twoway", headers = c("Poisson", "Gaussian"))

  gravity_subfe <- list()
  all_FEs <- c("Year", "Destination", "Origin")
  for (i in 0:3) {
    gravity_subfe[[i + 1]] <- fepois(Euros ~ log(dist_km), trade, fixef = all_FEs[0:i])
  }

  tab3 <- etable(gravity_subfe, cluster = ~ Origin + Destination)

  res_multi <- fepois(Euros ~ log(dist_km) | csw0(Year, Destination, Origin), trade)

  tab4 <- etable(res_multi, cluster = ~ Origin + Destination, tex = TRUE)

  myDict <- c("log(dist_km)" = "$\\ln (Distance)$", "(Intercept)" = "Constant")
  tab5 <- etable(res_multi,
    signifCode = c("a" = 0.01, "b" = 0.05),
    drop = "Const", dict = myDict, file = "Estimation Tables.tex",
    replace = TRUE, title = "First export -- normal Standard-errors"
  )
  tab6 <- etable(res_multi,
    cluster = ~Product, order = "Dist",
    dict = myDict, file = "Estimation Tables.tex",
    title = "Second export -- clustered standard-errors (on Product variable)"
  )

  gravity_pois_fes <- fixef(gravity_pois)
  s7 <- summary(gravity_pois_fes)

  gravity_pois_fes$Year
  plot(gravity_pois_fes)

  rm(gravity_pois, gravity_simple, gravity_pois_2, gravity_ols, gravity_negbin, gravity_subfe, res_multi)
  rm(s1, s2, s3, s4, s5, s6, s7)
  rm(tab1, tab2, tab3, tab4, tab5, tab6)
  rm(gravity_pois_fes)

  gc()
}

if (eval_quakes) {
  fit1 <- feols(depth ~ mag, quakes, "conley")
  fit2 <- feols(depth ~ mag, quakes, conley(200, distance = "spherical"))
  fit3 <- feols(depth ~ mag, quakes, vcov_conley(
    lat = "lat", lon = "long",
    cutoff = 200, distance = "spherical"
  ))

  rm(fit1, fit2, fit3)
  gc()
}

if (eval_did) {
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

  rm(est1, est2, est3, est4, est5, est6, est7, est8)
  rm(s1, s2, s3, s4)
  gc()
}

# SMALL SAMPLE CORRECTION ----

if (eval_smallsample) {
  est <- feols(y ~ x1 | id, base_did)
  est_up <- feols(y ~ x1 | id, base_did, ssc = ssc(fixef.K = "full"))
  est_down <- feols(y ~ x1 | id, base_did,
    ssc = ssc(adj = FALSE, cluster.adj = FALSE)
  )

  t1 <- etable(est, est_up, est_down)
  t2 <- etable(est, vcov = list(
    ~id, ~ id + ssc(fixef.K = "full"),
    ~ id + ssc(adj = FALSE, cluster.adj = FALSE)
  ))

  m1 <- feols(y ~ x1 | id, base_did, iid ~ ssc(adj = FALSE))
  m2 <- feols(y ~ x1 | id, base_did, hetero ~ ssc(adj = FALSE))

  rm(est, est_up, est_down, t1, t2, m1, m2)
  gc()
}

# OTHER VCOVS ----

if (eval_othervcovs) {
  est <- feols(y ~ x1 | id, base_did)
  summ <- summary(est, vcov = sandwich::vcovHC, type = "HC1")
  fit <- feols(y ~ x1 | id, base_did, vcov = function(x) sandwich::vcovHC(x, type = "HC1"))

  rm(est, summ, fit)
  gc()
}

# IVS ----

if (eval_ivs) {
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

  e1 <- etable(summary(est_iv_fe, stage = 1:2), fitstat = ~ . + ivfall +
    ivwaldall.p)

  rm(base, est_iv, est_iv_fe, fs, s1, e1)
  gc()
}

# INTERACTION TERMS ----

if (eval_interactions) {
  # base <- iris
  # names(base) <- c("y", paste0("x", 1:3), "fe1")
  # base$fe2 <- rep(letters[1:5], 30)
  # est_comb <- feols(y ~ x1 | fe1^fe2, base)
  # fes <- fixef(est_comb)[[1]]
  # est_vs <- feols(y ~ x1 | fe1[x2], base)
  # s_vs <- summary(fixef(est_vs))

  feols(Ozone ~ Solar.R + i(Month), airquality)

  # res_i1 <- feols(Ozone ~ Solar.R + i(Month), airquality)
  # res_i2 <- feols(Ozone ~ Solar.R + i(Month, ref = 8), airquality)
  # res_i3 <- feols(Ozone ~ Solar.R + i(Month, keep = 5:6), airquality)

  # t_i123 <- etable(res_i1, res_i2, res_i3,
  #   dict = c("6" = "June", "Month::5" = "May"),
  #   order = c("Int|May", "Mon")
  # )

  # est_did <- feols(y ~ x1 + i(period, treat, 5) | id + period, base_did)
  # iplot(est_did)

  # if (!require("ggplot2")) install.packages("ggplot2")
  # library(ggplot2)

  # g <- ggplot(
  #   aggregate(base_stagg[, c("year_treated", "treatment_effect_true")],
  #     by = list(
  #       year = base_stagg$year,
  #       group = to_integer(base_stagg$year_treated)
  #     ),
  #     mean
  #   ),
  #   aes(year, group, fill = year >= year_treated, alpha = treatment_effect_true)
  # ) +
  #   geom_tile(colour = "white", lwd = 1) +
  #   scale_fill_brewer("Treated?", palette = "Set1") +
  #   scale_alpha("Avg. treatment\neffect") +
  #   labs(x = "Year", y = "Group") +
  #   theme_minimal()

  # res_twfe <- feols(y ~ x1 + i(time_to_treatment, ref = c(-1, -1000)) |
  #   id + year, base_stagg)
  # res_sa20 <- feols(y ~ x1 + sunab(year_treated, year) | id + year, base_stagg)

  # iplot(list(res_twfe, res_sa20), sep = 0.5)

  # att_true <- tapply(
  #   base_stagg$treatment_effect_true,
  #   base_stagg$time_to_treatment, mean
  # )[-1]
  # points(-9:8, att_true, pch = 15, col = 4)

  # legend("topleft",
  #   col = c(1, 4, 2), pch = c(20, 15, 17),
  #   legend = c("TWFE", "Truth", "Sun & Abraham (2020)")
  # )

  # s_sa20 <- summary(res_sa20, agg = "att")
  # e_sa20 <- etable(res_sa20, agg = FALSE)

  # rm(
  #   base, est_comb, fes, est_vs, s_vs, res_i1, res_i2, res_i3, t_i123,
  #   est_did, g, res_twfe, res_sa20, att_true, s_sa20, e_sa20
  # )
  # gc()
}

# FORMULAS ----

if (eval_formulas) {
  base <- iris
  names(base) <- c("y", "x1", "x2", "x3", "species")

  setFixest_fml(..ctrl = ~ poly(x2, 2) + poly(x3, 2))

  xpd(y ~ x1 + ..ctrl)

  vars <- c("x2", "x2^2", "x3")
  for (i in 1:3) {
    print(xpd(y ~ x1 + ..ctrl, ..ctrl = vars[1:i]))
  }

  feols(y ~ x1 + ..ctrl, base)

  xpd(Armed.Forces ~ Population + regex("GNP|ployed"), data = longley)
}
