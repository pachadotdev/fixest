devtools::load_all()

# data(trade)
# gravity_pois = fepois(Euros ~ log(dist_km) | Origin + Destination + Product + Year, trade, debug = F)

# dataset <- read.csv("https://stats.idre.ucla.edu/stat/data/poisson_sim.csv")
# dataset$prog <- factor(dataset$prog,
#                        levels = 1:3,
#                        labels = c("General", "Academic", "Vocational")
# )
# saveRDS(dataset, "dev/poisson_sim.rds")

dataset <- readRDS("dev/poisson_sim.rds")

m1 <- glm(num_awards ~ 0 + prog + math, family = poisson, data = dataset)

m2 <- fepois(num_awards ~ math | prog, data = dataset)

m2$coefficients

# m2$env <- NULL
#
# print(m2)
#
# m2original <- readRDS("~/github/fixest2/dev/fit_fepois_original.rds")
#
# all.equal(m2, m2original)
#
# m2$wols$coefficients
# m2original$wols$coefficients
