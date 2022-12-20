X = cbind(1, mtcars$wt)
y = mtcars$mpg
w = 1
colnames(X) = c("(Intercept)","mpg")

cpp11::cpp_source("dev/matrix_r_cpp11.cpp")
Rcpp::sourceCpp("dev/matrix_r_rcpp.cpp")

a = matrix_r_rcpp(X)
b = matrix_r_cpp11(X)

a
b

all.equal(a,b)
