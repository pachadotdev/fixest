# devtools::clean_dll()
devtools::load_all()

Rcpp::sourceCpp("dev/correct_r_matrix.cpp")

# 2x2 ----

mtcars2 = mtcars[grepl("Merc 2", rownames(mtcars)),]
mtcars2 = mtcars2[1:3,]

X = cbind(1, mtcars2$wt)
y = mtcars2$mpg
w = 1
colnames(X) = c("(Intercept)","wt")

correct_0w = F
nthreads = 1
collin.tol = 10^(-10)

info_products = cpp_sparse_products(X, w, y, correct_0w, nthreads)
xwx = info_products$XtX
xwy = info_products$Xty
info_inv = cpp_cholesky(xwx, collin.tol, nthreads)

info_inv

# 3x3 ----

mtcars2 = mtcars[grepl("Merc 2", rownames(mtcars)),]
mtcars2 = mtcars2[1:3,]

X = cbind(1, mtcars2$wt, mtcars2$hp)
y = mtcars2$mpg
w = 1
colnames(X) = c("(Intercept)","wt", "hp")

correct_0w = F
nthreads = 1
collin.tol = 10^(-10)

info_products = cpp_sparse_products(X, w, y, correct_0w, nthreads)
xwx = info_products$XtX
xwy = info_products$Xty
info_inv = cpp_cholesky(xwx, collin.tol, nthreads)

info_inv
