#Matrix multiplication


pacman::p_load(igraph, intergraph, Matrix)


sparse <- as.matrix(as_adjacency_matrix(mtest, type = "both", sparse = T))
head(sparse)

adj.perc <- sparse/rowSums(sparse)

vec <- rep(T, 100)

adopters <- adj.perc %*% vec
sum(adopters)

testmat <- matrix(0,nrow = 5, ncol = 5)
testmat[1,2] <- 1
testmat[1,5] <- 1
testmat[2,1] <- 1
testmat[2,4] <- 3

vec <- c(T,F,F,T,F)



mtest <- make_lattice(c(10,10), nei = 2)

#mtest <- asNetwork(mtest)


adjmat <- as.matrix(as_adjacency_matrix(mtest, type = "both", sparse = T))
adj.perc <- adjmat/rowSums(adjmat)

n = 100
z = 10
meantau = 0.16
n_high = n/z

#Calculating the penalty to be subtracted from all nodes 
extra_influ = ((meantau * z) - 1) * n_high
penalty <- extra_influ / n / z

#Subtracting the penalty
adj.perc <- ifelse(adj.perc != 0, adj.perc - penalty, adj.perc)

#Sampling the high status nodes we want 
high_status_nodes <- sample(seq_len(n), size = n_high)

for (i in high_status_nodes) {
  node = i
  neighbors = get.neighborhood(asNetwork(mtest), node)
  adj.perc[neighbors,node] <- meantau
}

class(adj.perc[,10][1])

vec <- rep(F, 100)
vec[11] <- T

adopters <- adj.perc %*% vec

sum(adopters)
