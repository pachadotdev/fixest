library(tradepolicy)

ch1_application1 <- agtpa_applications %>%
    select(exporter, importer, pair_id, year, trade, dist, cntg, lang, clny) %>%
    filter(year %in% seq(1986, 2006, 4))

ch1_application1 <- ch1_application1 %>%
    mutate(
        log_trade = log(trade),
        log_dist = log(dist)
    )

ch1_application1 <- ch1_application1 %>%
    # Create Yit
    group_by(exporter, year) %>%
    mutate(
        y = sum(trade),
        log_y = log(y)
    ) %>%

    # Create Eit
    group_by(importer, year) %>%
    mutate(
        e = sum(trade),
        log_e = log(e)
    )

ch1_application1 <- ch1_application1 %>%
    # Replicate total_e
    group_by(exporter, year) %>%
    mutate(total_e = sum(e)) %>%
    group_by(year) %>%
    mutate(total_e = max(total_e)) %>%

    # Replicate rem_exp
    group_by(exporter, year) %>%
    mutate(
        remoteness_exp = sum(dist *  total_e / e),
        log_remoteness_exp = log(remoteness_exp)
    ) %>%

    # Replicate total_y
    group_by(importer, year) %>%
    mutate(total_y = sum(y)) %>%
    group_by(year) %>%
    mutate(total_y = max(total_y)) %>%

    # Replicate rem_imp
    group_by(importer, year) %>%
    mutate(
        remoteness_imp = sum(dist / (y / total_y)),
        log_remoteness_imp = log(remoteness_imp)
    )

ch1_application1 <- ch1_application1 %>%
    # This merges the columns exporter/importer with year
    mutate(
        exp_year = paste0(exporter, year),
        imp_year = paste0(importer, year)
    )

ch1_application1 <- ch1_application1 %>%
    filter(exporter != importer)

use_data(ch1_application1)

fixest2::fepois(trade ~ log_dist + cntg + lang + clny | exp_year + imp_year, data = ch1_application1)

# this output has wrong std error and t value
# > fixest2::fepois(trade ~ log_dist + cntg + lang + clny | exp_year + imp_year, data = ch1_application1)
# [1] "CHECK HERE"
# Poisson estimation, Dep. Var.: trade
# Observations: 28,152
# Fixed-effects: exp_year: 414,  imp_year: 414
# Standard-errors: Clustered (exp_year)
# Estimate Std. Error t value Pr(>|t|)
# log_dist -0.840927        Inf       0        1
# cntg      0.437443        Inf       0        1
# lang      0.247477        Inf       0        1
# clny     -0.222490        Inf       0        1
# ---
#     Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# Log-Likelihood: -2,194,536.1   Adj. Pseudo R2: 0.963391
# BIC:  4,397,587.8     Squared Cor.: 0.946097
