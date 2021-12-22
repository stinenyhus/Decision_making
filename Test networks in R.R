##Trying out the network tutorial

# Following best practices, we begin by setting up the workspace
rm(list = ls())
setwd("~/Uni/4. semester/Exam/SocCult")
library(pacman)
p_load(statnet, tidyverse, magrittr, ggnetwork, EpiModel) #Last two are dedicated diffusion
#install.packages("netdiffuseR") #This does not work


# Typically, the purpose of a diffusion simulation is to illustrate the effects
# of a new diffusion process, or a new type of network structure, on the
# ultimate results of a diffusion process.  As such diffusion simulations are 
# highly variable -- after all, the point for each one is to show something new!

# However, despite the variability, most diffusion simulations build on a
# common set of building blocks:
# - an underlying network structure,
# - a transmission process, and
# - a running timer.
# In essence, we are defining an algorithm, or a set of rules, for how people
# get "infected" with a social contagion, and looking at how those rules perform
# under differnt conditions.

# To simulate random mixing, we construct the following rules:
# 1. Set up a world with a fixed number of people and a small fraction of 
#    initial infections, arranged randomly as a network
# 2. Pick an edge at random
# 3. If the edge is discordant (a susceptible-infected connection), flip a coin
#    to determine whether the non-infected person gets infected
# 4. Repeat steps 1 and 2 until everyone is infected, or until a certain number
#    of steps has passed


###COMPLEX CONTAGIONS 

### Set parameters in advance ###
n <- 100
tau <- 2 / 6
max.time <- 50
sqrt(n)


library(pacman)
p_load(intergraph, igraph, ggnetwork )

# Construct a small world network
set.seed(919)
sw.net <- igraph::watts.strogatz.game(1, 50, 2, p = 0.005) %>%
  intergraph::asNetwork(.) 
plot(sw.net)

igraph::degree.distribution()

directed <- igraph::as.directed(sw.net)
directed <- asNetwork(directed)
plot(directed)
triads <- triad.census(directed)
triads
is_weighted(directed)

adj.mat <- sw.net[, ]
diag(adj.mat) <- 0
weight <- graph_from_adjacency_matrix(adj.mat, mode = "undirected", weighted = TRUE, diag = FALSE)
weight <- asNetwork(weight)
weight$mel[[100]]$atl$weight

plot(weight)

#Making Moore lattice without rewiring
moore <- make_lattice(c(100,100), nei =2) %>% asNetwork(.)
plot(moore)
moorenat <- moore[,]

#Making Moore lattice with rewiring
moore_re <- make_lattice(c(sqrt(n),sqrt(n)), nei =2) %>% rewire(each_edge(p = .005)) %>% asNetwork(.)
moore_re
plot(moore_re)

sqrt(10000)
490^2
sqrt(50)
7*7

#Random degree - how???? 

plot(sw.net)

#Maybe making network with power law neighbors???
degs <- sample(1:50, 50, replace=TRUE, prob=(1:50)^-1.5)
if (sum(degs) %% 2 != 0 ) {degs[1] <- degs[1] + 1}
g5 <- sample_degseq(degs, method="vl") #%>% asNetwork(.)
all(degree(g5) == degs)
g5 <- g5 %>% rewire(each_edge(p = 0)) %>% asNetwork(.)

plot(g5)

### Step 2: Choose a person & neighbors at random as an initial adopter ###
# Set up vector to indicate adoption
adopters <- rep(F, 50)
adopters

# Choose a person at random
initial.adopter <- sample(seq_len(n), size = 1)

# Get the list of people they're attached to
initial.neighbors <- get.neighborhood(sw.net, initial.adopter)
initial.neighbors


# Set them all as "adopters"
adopters[c(initial.adopter, initial.neighbors)] <- T
adopters

### Step 3: infect all people who have more than tau neighbors infected ###


# We can update everyone simultaneously using matrix multiplication
adj.mat <- sw.net[, ]
diag(adj.mat) <- 0  # set the diagonal to 0, b/c people don't weight themselves

adj.mat <- moore_re[, ]
diag(adj.mat) <- 0 

adj.mat <- g5[, ]
adj.mat
diag(adj.mat) <- 0 

# We could take the sum...
adj.mat %*% adopters

#... or the percentage
(adj.mat.rn <- adj.mat / rowSums(adj.mat))
246/99/4
0.25-0.0062
0.62*4*99 + 4

p <- adj.mat.rn %*% adopters
p

# And then we calculate the people who are above our threshold
(adj.mat.rn %*% adopters) >= tau

# Note that this actually allows people to abandon the innovation if enough of
# of their neighbors are not adopters.  For now we don't want that to happen, 
# so we'll only test people who are not yet adopters
ifelse(adopters, TRUE, ((adj.mat.rn %*% adopters) >= tau))

### Step 4: Repeat ###
# Again, we can take care of this by wrapping it in a loop
adopt <- vector(mode = "list", length = max.time)
adopt[[1]] <- adopters
adopt[[1]]

for (t in 2:max.time) {
  adopt[[t]] <- ifelse(adopters, TRUE, ((adj.mat.rn %*% adopt[[t - 1]]) >= tau))
}
id = 1:n
id

#Stochastic threshold 
p <- 5/10
M <- 10
L = 1/(1+exp((.5-p)*M))
L
activation <- rbinom(1,1,prob = L)
if activation == 1 then TRUE, else FALSE


if stochastic == TRUE 
elif stochastic == FALSE 
ifelse(adopters, TRUE, ((adj.mat.rn %*% adopt[[t - 1]]) >= tau))


# Note that again we get the characteristic S-shaped curve:
data.frame(
  t = 1:max.time,
  n.adopt = sapply(adopt, sum)
) %>%
  ggplot(aes(x = t, y = n.adopt)) + 
  geom_line()

set.seed(330)
g5.net.layout <- ggnetwork(g5) %>% 
  rename(id = vertex.names)

g5.net.layout.by.time <- adopt %>%
  lapply(FUN = as.data.frame) %>% 
  lapply(FUN = set_names, value = "adopter") %>% 
  lapply(FUN = mutate, id = 1:50) %>%
  lapply(FUN = right_join, y = sw.net.layout, by = "id") %>% 
  bind_rows(.id = "t") %>% 
  mutate(t = as.integer(t))

g5.net.layout.by.time %>% 
  filter(t > 30 & t < 40) %>% 
  ggplot(aes(xend = xend, yend = yend, x = x, y = y)) + 
  geom_edges(color = "lightgray") +
  geom_nodes(aes(color = adopter)) + 
  facet_wrap(~ t) + 
  theme_blank()

p_load(boot)

#Stochastic threshold 
p <- 5/10
M <- 10
L = 1/(1+exp((.5-p)*M))
L

rbinom(1,1,prob = L)


#Sample from binomial like this?
sample(c(0,1), 1, prob = prob)
#Or like this
rbinom(1,1, prob = prob)
