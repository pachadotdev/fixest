library(dplyr)
library(tidyr)

ch1_application3 <- agtpa_applications %>%
    filter(year %in% seq(1986, 2006, 4)) %>%
    mutate(
        exp_year = paste0(exporter, year),
        imp_year = paste0(importer, year),
        year = paste0("intl_border_", year),
        log_trade = log(trade),
        log_dist = log(dist),
        intl_brdr = ifelse(exporter == importer, pair_id, "inter"),
        intl_brdr_2 = ifelse(exporter == importer, 0, 1),
        pair_id_2 = ifelse(exporter == importer, "0-intra", pair_id)
    ) %>%
    spread(year, intl_brdr_2, fill = 0)

ch1_application3 <- ch1_application3 %>%
    group_by(pair_id) %>%
    mutate(sum_trade = sum(trade)) %>%
    ungroup()

form <- trade ~ 0 + rta + rta_lag4 + rta_lag8 + rta_lag12 +
    intl_border_1986 + intl_border_1990 + intl_border_1994 +
    intl_border_1998 + intl_border_2002 |
    exp_year + imp_year + pair_id_2

d <- filter(ch1_application3, sum_trade > 0)

out <- list()

for (i in 1:100) {
    print(i)
    out[[i]] <- fixest2::feols(form, data = d)$coefficients
}

for (i in 2:100) {
    x <- all.equal(
        out[[1]],
        out[[i]]
    )

    print(paste(i, x))
}

out[[84]]
