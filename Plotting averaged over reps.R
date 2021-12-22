#Running test ABMs of all conditions with fewer nodes and only one repetition

source("useful_function_soccult.r")

library(pacman)
p_load(tidyverse, network, igraph, intergraph, tidygraph, ggraph, ggplot2)


#Baseline
baseline <- read.csv("basetau_neu_rep100_tau25.csv")
baseline %>% group_by(round) %>% summarise(N_adopters = mean(adopters)) %>% filter(round < 10) %>% 
  ggplot(aes(round, N_adopters))+
  geom_line()


plot_standard(test, title = "Baseline")

#Baseline with 
base_10000_r01_neu_rep10 <- contagion_sim(rep = 10, rounds = 50, tau = .25, n = 10000, nei = 2, p = 0.1)
write.csv(base_10000_r01_neu_rep10, "base_10000_r01_neu_rep10.csv")

#Baseline with 22500 - startet 17:03
base_22500_r01_neu_rep10 <- contagion_sim(rep = 1, rounds = 50, tau = .25, n = 22500, nei = 2, p = 0.1)
base_22500_r01_neu_rep10 %>% 
  group_by(round) %>% 
  summarise(sum_adopt = sum(adopters)) %>%
  ggplot(aes(x = round, y = sum_adopt)) + 
  geom_line()

#Baseline with 22500 - startet 17:03
base_22500_r01_neu_rep10_100rounds <- contagion_sim(rep = 1, rounds = 100, tau = .25, n = 22500, nei = 2, p = 0.1)
base_22500_r01_neu_rep10_100rounds %>% 
  group_by(round) %>% 
  summarise(sum_adopt = sum(adopters)) %>%
  ggplot(aes(x = round, y = sum_adopt)) + 
  geom_line()



#high status
base_100_r01_neu_highstatus <- contagion_sim(high_node = T, rep = 1, rounds = 50, tau = .25, n = 100, nei = 2, p = 0.1)
write.csv(base_100_r01_neu_highstatus, "base_100_r01_neu_highstatus.csv")

base_100_r01_neu_highstatus %>% 
  group_by(round) %>% 
  summarise(sum_adopt = sum(adopters)) %>%
  ggplot(aes(x = round, y = sum_adopt)) + 
  geom_line()

#stochastic
base_100_r01_neu_stochastic <- contagion_sim(rep = 1, rounds = 100, tau = .25, n = 100, nei = 2, p = 0.1, stochastic = T)
write.csv(base_100_r01_neu_stochastic, "base_100_r01_neu_stochastic")

base_100_r01_neu_stochastic %>% 
  group_by(round) %>% 
  summarise(sum_adopt = sum(adopters)) %>%
  ggplot(aes(x = round, y = sum_adopt)) + 
  geom_line()

#scale free
base_100_scalefree <- contagion_sim(rep = 1, rounds = 100, tau = .25, net_type = "scale_free", degrees = 50, gamma = 2.3, n = 100, nei = 2, p = 0)
write.csv(base_100_scalefree, "base_100_scalefree")

base_100_scalefree %>% 
  group_by(round) %>% 
  summarise(sum_adopt = sum(adopters)) %>%
  ggplot(aes(x = round, y = sum_adopt)) + 
  geom_line()

#varying thresholds
randomtau_100_r01_neu <- contagion_sim(tau_type = "random_tau", rep = 1, rounds = 100, tau = .25, n = 100, nei = 2, p = 0.1) 
write.csv(randomtau_100_r01_neu, "randomtau_100_r01_neu")

randomtau_100_r01_neu %>% 
  group_by(round) %>% 
  summarise(sum_adopt = sum(adopters)) %>%
  ggplot(aes(x = round, y = sum_adopt)) + 
  geom_line()


#stochastic x varyig thresholds
randomtau_100_r01_neu_stochastic <- contagion_sim(tau_type = "random_tau", rep = 1, rounds = 100, tau = .25, n = 100, nei = 2, p = 0.1, stochastic = T) 
write.csv(randomtau_100_r01_neu_stochastic, "randomtau_100_r01_neu_stochastic")

randomtau_100_r01_neu_stochastic %>% 
  group_by(round) %>% 
  summarise(sum_adopt = sum(adopters)) %>%
  ggplot(aes(x = round, y = sum_adopt)) + 
  geom_line()


#scalefree x varying thresholds
randomtau_100_r01_scalefree <- contagion_sim(tau_type = "scale_free", degrees = 50, gamma = 2.3, rep = 1, rounds = 100, tau = .25, n = 100, nei = 2, p = 0.1) 
write.csv(randomtau_100_r01_scalefree, "randomtau_100_r01_scalefree")

#100 reps baseline
base_100_r01_neu_100rep <- contagion_sim(rep = 100, rounds = 100, tau = .25, n = 100, nei = 2, p = 0.1)

#high status 22500 nodes- startet at 15:25
base_22500_r01_neu_highstatus <- contagion_sim(high_node = T, rep = 1, rounds = 100, tau = .25, n = 22500, nei = 2, p = 0.1)


#Averaging over rounds 
base_100_r01_neu_100rep %>% group_by(network,round) %>% 
  summarise(sum_adopt = sum(adopters)) %>% 
  group_by(round) %>% summarise(mean = mean(sum_adopt)) %>%
  ggplot(aes(x = round, y = mean)) +
  geom_line()


#Loading files
test <- read.csv("basetau_neu_rep1_tau25.csv")

test %>% 
  group_by(round) %>% summarise(mean = mean(adopters)) %>%
  ggplot(aes(x = round, y = mean)) +
  geom_line()

#Plotting from different dataframes
ggplot()+
  geom_line(data=simulation1_summed, aes(round, sumadopt, color = "simulation 1"), size = 1.2)+
  geom_line(data=simulation3_summed, aes(round, sumadopt, color = "simulation 2"), size = 1.2)+
  scale_colour_pander() +
  scale_fill_pander() +
  theme_pander()
