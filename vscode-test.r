# devtools::clean_dll()
devtools::load_all()

# Rcpp::sourceCpp("dev/correct_r_matrix.cpp")

test_chol_2x2 <- T
test_chol_3x3 <- T

test_simple_ols <- T

# Test Choslesky 2x2 ----

if (test_chol_2x2) {
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
}

# 3x3 ----

if (test_chol_3x3) {
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
}

# simple LM ----

if (test_simple_ols) {
    fit1 <- feols(mpg ~ wt, data = mtcars)
    fit2 <- lm(mpg ~ wt, data = mtcars)

    all.equal(fit1$coefficients, fit2$coefficients)
    all.equal(fit1$residuals, unname(fit2$residuals))
}

# weighted LM ----

set.seed(0)
base <- iris
names(base) <- c("y", "x1", "x2", "x3", "species")
base$offset_value <- unclass(base$species) - 0.95
base$fe_2 <- rep(1:5, 30)

# works
feols(y ~ x1, data = base)
summary(lm(y ~ x1, data = base))

# fails
feols(y ~ x1 | fe_2, data = base)
summary(lm(y ~ 0 + x1 + as.factor(fe_2), data = base))
summary(lfe::felm(y ~ x1 | fe_2, data = base))

