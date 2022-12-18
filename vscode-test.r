library(fixest)
check <- feols(mpg ~ wt, data = mtcars)
