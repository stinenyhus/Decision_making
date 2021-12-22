#Running the models
setwd("/work/Exam")
source("useful_function_soccult.r")
install.packages("pacman")
library(pacman)
p_load(tidyverse, network, igraph, intergraph, tidygraph, ggraph, ggplot2, ggthemes)

test_small_net <- contagion_sim(tau_type = "random_tau",
                                      tau = 0.33, 
                                      n = 100,
                                      nei = 2,
                                      rep = 5, 
                                      rounds = 100)

# This test was 10 connected nodes, 50 connectedness, tau=0.25, n = 22500, 5 rep of 100 rounds
test_connected_nodes <- read.csv("test_connected_nodes_2112.csv")

test_small_net %>%
  ggplot(aes(round, adopters))+
  geom_line()+
  facet_wrap(~network)
head(basetau_neu_rep100_tau25)

test_summed <- sum_data(test_connected_nodes, name="test", tau="0.25")
plot_standard_by_name(test_summed, "Test run")

net <- generate_neumann(22500, 0, 2)
adj <- as.matrix(as_adjacency_matrix(net, type = "both", sparse = T))
diag(adj) <- 0 
adj = create_connected_nodes(adj, 10, 50)
graph = graph_from_adjacency_matrix(adj)
net = asNetwork(graph)
asNetwork(adj)

df = data.frame(node=rep(0,22500), neighbors=rep(0,22500), n_neighbors = rep(0,22500))
for (i in 1:22500){
  neighbors = get.neighborhood(net, i)
  df$node[i] = i
  df$neighbors[i] = toString(neighbors)
  df$n_neighbors[i] = length(neighbors)
}
df %>% subset(n_neighbors > 40)

df %>% group_by(n_neighbors) %>% summarize(n=n())

plot.igraph(graph)
####STINES MODELS####

#Baseline - startet 10:36 finished at 16:20
basetau_neu_rep100_tau25 <- contagion_sim(rep = 100, rounds = 100, tau = .25, n = 22500, nei = 2, p = 0.1)
write.csv(basetau_neu_rep100_tau25, "basetau_neu_rep100_tau25.csv")

#Startet at 18:17
basetau_neu_rep100_tau30 <- contagion_sim(rep = 100, rounds = 100, tau = .30, n = 22500, nei = 2, p = 0.1)
write.csv(basetau_neu_rep100_tau30, "basetau_neu_rep100_tau30.csv")

#Tau .40 baseline
basetau_neu_rep100_tau40 <- contagion_sim(rep = 100, rounds = 100, tau = .40, n = 22500, nei = 2, p = 0.1)
write.csv(basetau_neu_rep100_tau40, "basetau_neu_rep100_tau40.csv")

#Baseline - startet 10:09 finished at 16:10
basetau_neu_rep3_tau33_re1 <- contagion_sim(rep = 3, rounds = 100, tau = 4/12, n = 22500, nei = 2, p = 0.1)
basetau_neu_rep3_tau33_re01 <- contagion_sim(rep = 3, rounds = 100, tau = 4/12, n = 22500, nei = 2, p = 0.01)
basetau_neu_rep3_tau33_re05 <- contagion_sim(rep = 3, rounds = 100, tau = 4/12, n = 22500, nei = 2, p = 0.05)
basetau_neu_rep3_tau33_re03 <- contagion_sim(rep = 3, rounds = 100, tau = 4/12, n = 22500, nei = 2, p = 0.03)

#With highstatus nodes
basetau_neu_rep100_highstatus_tau25 <- contagion_sim(rep = 100, rounds = 100, tau = .25, n = 22500, nei = 2, p = 0.1, high_node = T)
write.csv(basetau_neu_rep100_highstatus, "basetau_neu_rep100_highstatus.csv")


####ASTRIDS MODELS####

#Stochastic
basetau_neu_rep100_stochastic <- contagion_sim(rep = 100, rounds = 100, tau = .25, n = 22500, nei = 2, p = 0.1, stochastic = T)
write.csv(basetau_neu_rep100_stochastic, "basetau_neu_rep100_stochastic.csv")

#Stochastic with other taus - startet ca 10:10 - 20 reps 12:20, 33 reps at 13:00, 50 at 14:00, 70 kl 15:45
basetau_neu_rep100_stochastic_tau36 <- contagion_sim(rep = 100, rounds = 100, tau = .36, n = 22500, nei = 2, p = 0.1, stochastic = T)
write.csv(basetau_neu_rep100_stochastic_tau36, "basetau_neu_rep100_stochastic_tau36.csv")
basetau_neu_rep100_stochastic_tau36_summed <- sum_data(basetau_neu_rep100_stochastic_tau36)
plot_standard(basetau_neu_rep100_stochastic_tau36_summed)

#Stochastic tau .4
basetau_neu_rep100_stochastic_tau40 <- contagion_sim(rep = 100, rounds = 100, tau = .4, n = 22500, nei = 2, p = 0.1, stochastic = T)
write.csv(basetau_neu_rep100_stochastic_tau40, "basetau_neu_rep100_stochastic_tau40.csv")
basetau_neu_rep100_stochastic_tau40_summed <- sum_data(basetau_neu_rep100_stochastic_tau40)
plot_standard(basetau_neu_rep100_stochastic_tau40_summed)

#Scalefree with minimally complex
basetau_scalefree_rep100_tau17 <- contagion_sim(net_type = "scale_free", rep = 100, rounds = 100, tau = 2/12, n = 22500, nei = 2, p = 0, degrees = 75, gamma = 2.3)
write.csv(basetau_scalefree_rep100_tau17, "basetau_scalefree_rep100_tau17.csv")
basetau_scalefree_rep100_tau17_summed <- sum_data(basetau_scalefree_rep100_tau17)
plot_standard(basetau_scalefree_rep100_tau17_summed)

####MARIES MODELS####

#Varying thresholds
randomtau_neu_rep1 <- contagion_sim(tau_type = "random_tau", rep = 1, rounds = 50, tau = .25, n = 22500, nei = 2, p = 0.1)
write.csv(randomtau_neu_rep100, "randomtau_neu_rep100.csv")

basetau_neu_rep100_tau30 <- contagion_sim(rep = 100, rounds = 100, tau = .30, n = 22500, nei = 2, p = 0.1)
write.csv(basetau_neu_rep100_tau30, "basetau_neu_rep100_tau30.csv")

#High status + stochastic
basetau_neu_rep1_highstatus_stochastic_22500 <- contagion_sim(rep = 1, rounds = 50, tau = .25, n = 22500, nei = 2, p = 0.1, stochastic = T, high_node = T)
write.csv(basetau_neu_rep100_highstatus_stochastic, "basetau_neu_rep100_highstatus_stochastic.csv")



####Extra models####

#Scalefree + varying thresholds
randomtau_scalefree_rep100_tau25 <- contagion_sim(net_type = "scale_free",tau_type = "random_tau", rep = 100, rounds = 100, tau = .25, n = 22500, nei = 2, p = 0, degrees = 75, gamma = 2.3)
write.csv(randomtau_scalefree_rep100_tau25, "randomtau_scalefree_rep100_tau25.csv")


#High status + stochastic + random tau
randomtau_neu_rep100_stochastic <- contagion_sim(tau_type = "random_tau", rep = 100, rounds = 100, tau = .25, n = 22500, nei = 2, p = 0.1, stochastic = T)
write.csv(randomtau_neu_rep100_highstatus_stochastic, "randomtau_neu_rep100_highstatus_stochastic.csv")

#Stochastic + varying tau (mean = .36) startet 9:30
randomtau_neu_rep100_tau36_stochastic <- contagion_sim(tau_type = "random_tau", rep = 100, rounds = 100, tau = .36, n = 22500, nei = 2, p = 0.1, stochastic = T)
write.csv(randomtau_neu_rep100_tau36_stochastic, "randomtau_neu_rep100_tau36_stochastic.csv")

#Random + stochastic tau, startet at 07:48 - 20 reps 10:05, 65 kl 14:35, 72 kl 15:48, 80 kl 16:36
randomtau_neu_rep100_tau35_stochastic <- contagion_sim(tau_type = "random_tau", rep = 100, rounds = 100, tau = .35, n = 22500, nei = 2, p = 0.1, stochastic = T)
write.csv(randomtau_neu_rep100_tau35_stochastic, "randomtau_neu_rep100_tau35_stochastic.csv")

#Stochastic + varying tau (mean = .4)
randomtau_neu_rep100_tau40_stochastic <- contagion_sim(tau_type = "random_tau", rep = 100, rounds = 100, tau = .40, n = 22500, nei = 2, p = 0.1, stochastic = T)
write.csv(randomtau_neu_rep100_tau40_stochastic, "randomtau_neu_rep100_tau40_stochastic.csv")
