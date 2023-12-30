input <- readRDS("dev/input.rds")

# explicitly convert to the right data types

Q <- as.integer(input$Q)

N <- as.integer(input$N)

S <- as.double(input$S)

dumMat <- input$dumMat
storage.mode(dumMat) <- "integer"

nbCluster <- as.integer(input$nbCluster)

orderCluster <- input$orderCluster
storage.mode(orderCluster) <- "integer"

# check the data types
print(paste("Q", class(Q)))
print(paste("N", class(N)))
print(paste("S", class(S)))
print(paste("dumMat", class(dumMat)))
print(paste("nbCluster", class(nbCluster)))
print(paste("orderCluster", class(orderCluster)))

# check dimensions
print(paste("Q", length(Q)))
print(paste("N", length(N)))
print(paste("S", length(S)))
print(paste("dumMat", paste("rows", dim(dumMat)[1], "cols", dim(dumMat)[2])))
print(paste("nbCluster", length(nbCluster)))
print(paste("orderCluster", paste("rows", dim(orderCluster)[1], "cols", dim(orderCluster)[2])))

# cpp11::cpp_source("dev/00_cpp11_get_fe_gnl.cpp")
cpp11::cpp_source("dev/00_cpp11_get_fe_gnl_2.cpp")

cpp11_get_fe_gnl(Q, N, S, dumMat, nbCluster, orderCluster)
