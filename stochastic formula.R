library(pacman)
p_load(tidyverse, network, igraph, intergraph, tidygraph, ggraph, ggplot2, ggthemes)

####Testing different stochastic functions - a mess ####
M = 10

#With old function 
tau = 0.5
p = seq(0, 1, by = 0.01)
M=10
q = 0.1-tau + 0.5
L = 1/(1+exp((0.5-q)*M))
plot(L)

tau = 0.25
q = 0-tau + 0.5
(1/(1+exp((0.5-q)*M)))/70

22.5/22500

#With new function
x = seq(0,1, by = 0.01) #Proportion of neighborhood infected
mu = 0.5#Changing this term changes where L = 0.5 
#But it also changes the curve 
#But is this just what needs changing? i.e. should this follow tau?
s = 10 #a term partly defining the curve

cdf=1/(1+exp(-(x-mu)/s))
plot(cdf)

#If x = 0 and mu = 0.25, then L = 0.07 which does not solve our problem
22500*0.02

M = 10
x= 0.1
L = 1/(1+exp((.5-x)*M))
plot(L)
rbinom(1,1, prob = L)



ggplot(sto, aes(p,L))+
  geom_line()+
  xlab("p = proportion of activated neighborhood")+
  ylab("L = probability of adopting")



d <- filter(simulation1, adopters >= 22500*0.25)
dd <- d %>% group_by(network) %>% summarise(minround = min(round))
round(mean(dd$minround),0)
sd(dd$minround)
d %>% group_by(round) %>% summarise(sumadopt = mean(adopters))

test1 <- rbinom(100,1, prob = 0.5)

rbinom(1,prob = 0.5)


source("useful_function_soccult.r")

###Testing different taus in stochastic formula 
basetau_neu_n2500_t25_stochastic <- contagion_sim(rep = 1, rounds = 100, tau = 3/12, n = 2500, nei = 2, p = 0.1, stochastic = T)
basetau_neu_n2500_t25_stochastic_summed <- sum_data(basetau_neu_n2500_t25_stochastic)
plot_standard(basetau_neu_n2500_t25_stochastic_summed)

####Interesting range of taus####

basetau_neu_n2500_t33_stochastic <- contagion_sim(rep = 10, rounds = 100, tau = .33, n = 2500, nei = 2, p = 0.1, stochastic = T)
basetau_neu_n2500_t33_stochastic_summed <- sum_data(basetau_neu_n2500_t33_stochastic, name = "t33")

basetau_neu_n2500_t33_stochastic_no_random <- contagion_sim(rep = 10, rounds = 100, tau = .33, n = 2500, nei = 2, p = 0.1, stochastic = T)
basetau_neu_n2500_t33_stochastic_no_random_summed <- sum_data(basetau_neu_n2500_t33_stochastic_no_random, name = "t33_no_random")
#
basetau_neu_n2500_t35_stochastic <- contagion_sim(rep = 10, rounds = 100, tau = .35, n = 2500, nei = 2, p = 0.1, stochastic = T)
basetau_neu_n2500_t35_stochastic_summed <- sum_data(basetau_neu_n2500_t35_stochastic, name = "t35")

basetau_neu_n2500_t35_stochastic_no_random <- contagion_sim(rep = 10, rounds = 100, tau = .35, n = 2500, nei = 2, p = 0.1, stochastic = T)
basetau_neu_n2500_t35_stochastic_no_random_summed <- sum_data(basetau_neu_n2500_t35_stochastic_no_random, name = "t35_no_random")
#
basetau_neu_n2500_t36_stochastic <- contagion_sim(rep = 10, rounds = 100, tau = .36, n = 2500, nei = 2, p = 0.1, stochastic = T)
basetau_neu_n2500_t36_stochastic_summed <- sum_data(basetau_neu_n2500_t36_stochastic, name = "t36")

basetau_neu_n2500_t36_stochastic_no_random <- contagion_sim(rep = 10, rounds = 100, tau = .36, n = 2500, nei = 2, p = 0.1, stochastic = T)
basetau_neu_n2500_t36_stochastic_no_random_summed <- sum_data(basetau_neu_n2500_t36_stochastic_no_random, name = "t36_no_random")
#

basetau_neu_n2500_t38_stochastic <- contagion_sim(rep = 10, rounds = 100, tau = .38, n = 2500, nei = 2, p = 0.1, stochastic = T)
basetau_neu_n2500_t38_stochastic_summed <- sum_data(basetau_neu_n2500_t38_stochastic, name = "t38")

basetau_neu_n2500_t38_stochastic_no_random <- contagion_sim(rep = 10, rounds = 100, tau = .38, n = 2500, nei = 2, p = 0.1, stochastic = T)
basetau_neu_n2500_t38_stochastic_no_random_summed <- sum_data(basetau_neu_n2500_t38_stochastic_no_random, name = "t38_no_random")
#

basetau_neu_n2500_t40_stochastic <- contagion_sim(rep = 10, rounds = 100, tau = .40, n = 2500, nei = 2, p = 0.1, stochastic = T)
basetau_neu_n2500_t40_stochastic_summed <- sum_data(basetau_neu_n2500_t40_stochastic, name = "t40")

basetau_neu_n2500_t40_stochastic_no_random <- contagion_sim(rep = 10, rounds = 100, tau = .40, n = 2500, nei = 2, p = 0.1, stochastic = T)
basetau_neu_n2500_t40_stochastic_no_random_summed <- sum_data(basetau_neu_n2500_t40_stochastic_no_random, name = "t40_no_random")


t33_to_t40 <- rbind(basetau_neu_n2500_t33_stochastic_summed, basetau_neu_n2500_t35_stochastic_summed, basetau_neu_n2500_t36_stochastic_summed,
                    basetau_neu_n2500_t38_stochastic_summed, basetau_neu_n2500_t40_stochastic_summed)

t33_to_t40_no_random <- rbind(basetau_neu_n2500_t33_stochastic_no_random_summed, basetau_neu_n2500_t35_stochastic_no_random_summed, 
                              basetau_neu_n2500_t36_stochastic_no_random_summed,
                              basetau_neu_n2500_t38_stochastic_no_random_summed, basetau_neu_n2500_t40_stochastic_no_random_summed)

plot_standard(t33_to_t40)
plot_standard(t33_to_t40_no_random)


####Too far out taus ####

basetau_neu_n2500_t42_stochastic <- contagion_sim(rep = 10, rounds = 100, tau = .42, n = 2500, nei = 2, p = 0.1, stochastic = T)
basetau_neu_n2500_t42_stochastic_summed <- sum_data(basetau_neu_n2500_t42_stochastic)
plot_standard(basetau_neu_n2500_t42_stochastic_summed)

basetau_neu_n2500_t44_stochastic <- contagion_sim(rep = 10, rounds = 100, tau = .44, n = 2500, nei = 2, p = 0.1, stochastic = T)
basetau_neu_n2500_t44_stochastic_summed <- sum_data(basetau_neu_n2500_t44_stochastic)
plot_standard(basetau_neu_n2500_t44_stochastic_summed)

basetau_neu_n2500_t46_stochastic <- contagion_sim(rep = 10, rounds = 100, tau = .46, n = 2500, nei = 2, p = 0.1, stochastic = T)
basetau_neu_n2500_t46_stochastic_summed <- sum_data(basetau_neu_n2500_t46_stochastic)
plot_standard(basetau_neu_n2500_t46_stochastic_summed)

basetau_neu_n2500_t50_stochastic <- contagion_sim(rep = 10, rounds = 100, tau = .50, n = 2500, nei = 2, p = 0.1, stochastic = T)
basetau_neu_n2500_t50_stochastic_summed <- sum_data(basetau_neu_n2500_t50_stochastic)
plot_standard(basetau_neu_n2500_t50_stochastic_summed)
