devtools::load_all()

X = cbind(1, mtcars$wt)
y = mtcars$mpg
w = 1
colnames(X) = c("(Intercept)","mpg")

correct_0w = F
nthreads = 1
collin.tol = 10^(-5)

info_products = cpp_sparse_products(X, w, y, correct_0w, nthreads)
xwx = info_products$XtX
xwy = info_products$Xty
info_inv = cpp_cholesky(xwx, collin.tol, nthreads)

info_inv

solve(xwx)

feols(
    mpg ~ wt,
    data = mtcars
)

stats::lm(
    mpg ~ wt,
    data = mtcars
)
