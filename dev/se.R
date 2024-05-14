library(fixest2)

lhs = 'y1'
rhs = 'x3'
all_rhs = c('', 'x2', 'x3')
n_rhs = 1
k = 1

base = readRDS("dev/base.rds")

# try to create a memory leak by repetition
for (i in 1:3) {
  est_multi <- feols(c(y1, y2) ~ x1 + csw0(x2, x3) + x4 | species + fe2, base, fsplit = ~species)
}
# ok

# try to create a memory leak by repetition
for (i in 1:3) {
  mat <- coeftable(object = est_multi[[1]],
      vcov = NULL, ssc = NULL, cluster = "fe3",
      keep = NULL, drop = NULL, order = NULL)

  mat[, 2]
}
# problem
