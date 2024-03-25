# if (!require("devtools")) install.packages("devtools")
# devtools::load_all()

library(fixest2)

data(base_did)

# est_panel <- feols(y ~ x1, base_did, panel.id = ~ id + period)
# summary(est_panel, "newey_west")

base <- iris
names(base) <- c("y", "x1", "x_endo_1", "x_inst_1", "fe")
set.seed(2)
base$x_inst_2 <- 0.2 * base$y + 0.2 * base$x_endo_1 + rnorm(150, sd = 0.5)
base$x_endo_2 <- 0.2 * base$y - 0.2 * base$x_inst_1 + rnorm(150, sd = 0.5)

est_iv <- feols(y ~ x1 | x_endo_1 + x_endo_2 ~ x_inst_1 + x_inst_2, base)

# fixest 2
# > est_iv <- feols(y ~ x1 | x_endo_1 + x_endo_2 ~ x_inst_1 + x_inst_2, base)
# Error: Invalid input type, expected 'double' actual 'NULL'
# > est_iv
# Error: object 'est_iv' not found
# > est_iv <- feols(y ~ x1 | x_endo_1 + x_endo_2 ~ x_inst_1 + x_inst_2, base)
# Error: in feols(env = current_env, xwx = ZXtZX, xwy = ZXtu[...:
# All variables, '(Intercept)', 'x_inst_1' and 2 others, are virtually
# constant and equal to 0. Without doubt, your model is misspecified.
# > est_iv <- feols(y ~ x1 | x_endo_1 + x_endo_2 ~ x_inst_1 + x_inst_2, base)
# > est_iv
# TSLS estimation - Dep. Var.: y
#                   Endo.    : x_endo_1, x_endo_2
#                   Instr.   : x_inst_1, x_inst_2
# Second stage: Dep. Var.: y
# Observations: 150
# Standard-errors: NA (not-available)
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept)       NaN        NaN     NaN       NA
# fit_x_endo_1      NaN        NaN     NaN       NA
# fit_x_endo_2      NaN        NaN     NaN       NA
# x1                NaN        NaN     NaN       NA
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# F-test (1st stage), x_endo_1: stat = NA, p = NA, on 2 and 146 DoF.
# F-test (1st stage), x_endo_2: stat = NA, p = NA, on 2 and 146 DoF.
#                   Wu-Hausman: stat = NA, p = NA, on 2 and 144 DoF.

# est_iv <- fixest::feols(y ~ x1 | x_endo_1 + x_endo_2 ~ x_inst_1 + x_inst_2, base)
# > est_iv
# TSLS estimation - Dep. Var.: y
#                   Endo.    : x_endo_1, x_endo_2
#                   Instr.   : x_inst_1, x_inst_2
# Second stage: Dep. Var.: y
# Observations: 150
# Standard-errors: IID
#              Estimate Std. Error  t value   Pr(>|t|)
# (Intercept)  1.831380   0.411435  4.45121 1.6844e-05 ***
# fit_x_endo_1 0.444982   0.022086 20.14744  < 2.2e-16 ***
# fit_x_endo_2 0.639916   0.307376  2.08186 3.9100e-02 *
# x1           0.565095   0.084715  6.67051 4.9180e-10 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# RMSE: 0.398842   Adj. R2: 0.761653
# F-test (1st stage), x_endo_1: stat = 903.2    , p < 2.2e-16 , on 2 and 146 DoF.
# F-test (1st stage), x_endo_2: stat =   3.25828, p = 0.041268, on 2 and 146 DoF.
#                   Wu-Hausman: stat =   6.79183, p = 0.001518, on 2 and 144 DoF.
