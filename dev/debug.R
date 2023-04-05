load("~/github/fixest2/data/ch1_application1.rda")
fixest2::fepois(trade ~ log_dist + cntg + lang + clny | exp_year + imp_year, data = ch1_application1)
