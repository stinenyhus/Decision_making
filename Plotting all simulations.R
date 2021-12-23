



baseline_50 <- read.csv("data/NoHigh_tau_0.33_nSeeds_50/simulation_results.csv")
baseline_100 <- read.csv("data/NoHigh_tau_0.33_nSeeds_100/simulation_results.csv")
baseline_150 <- read.csv("data/NoHigh_tau_0.33_nSeeds_150/simulation_results.csv")

# High status, status = 0.5% of tau
status05_50 <-read.csv("data/nHigh_50_highStatus_TRUE_HighNodeTauPerc_0.5/simulation_results.csv")
status05_100 <-read.csv("data/nHigh_100_highStatus_TRUE_HighNodeTauPerc_0.5/simulation_results.csv")
status05_150 <-read.csv("data/nHigh_150_highStatus_TRUE_HighNodeTauPerc_0.5/simulation_results.csv")

# Summing
base_50_sum <- sum_data(baseline_50, name= "Baseline network with 50 initial seeds")
base_100_sum <- sum_data(baseline_100, name= "Baseline network with 100 initial seeds")
base_150_sum <- sum_data(baseline_150, name= "Baseline network with 150 initial seeds")

stat05_50_sum <- sum_data(status05_50, name = "0.5% of tau, 50 initial seed")
stat05_100_sum <- sum_data(status05_100, name = "0.5% of tau, 100 initial seed")
stat05_150_sum <- sum_data(status05_150, name = "0.5% of tau, 150 initial seed")

# Combining dfs
baseline <- rbind(base_50_sum, base_100_sum, base_150_sum)
status05 <- rbind(stat05_50_sum, stat05_100_sum, stat05_150_sum)

# Plotting
plot_standard_by_name(baseline, title = "Spread of contagion for baseline networks")
plot_standard_by_name(status05, title = "Spred of contagion for network where high status=0.5*tau")


ggplot(base_50_sum, aes(round, sumadopt))+ # color = str_wrap(tau,20)
  geom_line(size = 1.2)+
  theme_minimal()+
  theme(text = element_text(size = 15, family = "serif"), legend.key.height = unit(1,"cm"))+
  scale_color_brewer(palette = "PuOr")+
  labs(title = "Baseline network - 50 seeds")



calculate_point_estimates <- function(data, data_summed){
  q25=filter(data, adopters >= 22500*0.25) %>% group_by(network) %>% summarise(minround = min(round))
  q25_mean = round(mean(q25$minround),0)
  q25_sd   = round(sd(q25$minround),2)
  q50=filter(data, adopters >= 22500*0.50) %>% group_by(network) %>% summarise(minround = min(round))
  q50_mean = round(mean(q50$minround),0)
  q50_sd   = round(sd(q50$minround),2)
  q75=filter(data, adopters >= 22500*0.75) %>% group_by(network) %>% summarise(minround = min(round))
  q75_mean = round(mean(q75$minround),0)
  q75_sd   = round(sd(q75$minround),2)
  q99=filter(data, adopters >= 22500*0.99) %>% group_by(network) %>% summarise(minround = min(round))
  q99_mean = round(mean(q99$minround),0)
  q99_sd   = round(sd(q99$minround),2)
  all <- data.frame(
    name = data_summed$name[1],
    mean_25 = q25_mean,
    sd_25   = q25_sd,
    mean_50 = q50_mean,
    sd_50   = q50_sd,
    mean_75 = q75_mean,
    sd_75   = q75_sd,
    mean_99 = q99_mean,
    sd_99   = q99_sd
  )
  return(all)
}


































##### SOCCULT #####
#
####Simulation 0 - baseline####
#Loading file 
simulation0_all_sum <- read.csv("simulation0_all_sum.csv")

simulation0_all_sum$tau <- ifelse(simulation0_all_sum$tau == "0.25", ".25",simulation0_all_sum$tau)
simulation0_all_sum$tau <- ifelse(simulation0_all_sum$tau == "0.3", ".30",simulation0_all_sum$tau)
simulation0_all_sum$tau <- ifelse(simulation0_all_sum$tau == "0.35", ".35",simulation0_all_sum$tau)
simulation0_all_sum$tau <- ifelse(simulation0_all_sum$tau == "0.4", ".40",simulation0_all_sum$tau)

#Plotting
simulation0_plot <- plot_standard(simulation0_all_sum, 
                                  title = "Simulation 0 - baseline with varying tau")

#
####Simulation 1 - heterogeneity of thresholds####
#Loading files
simulation1_all_sum <- read.csv("simulation1_all_sum.csv")
simulation1_all_sum$tau <- ifelse(simulation1_all_sum$tau == "0.25", ".25",simulation1_all_sum$tau)
simulation1_all_sum$tau <- ifelse(simulation1_all_sum$tau == "0.3", ".30",simulation1_all_sum$tau)
simulation1_all_sum$tau <- ifelse(simulation1_all_sum$tau == "0.35", ".35",simulation1_all_sum$tau)
simulation1_all_sum$tau <- ifelse(simulation1_all_sum$tau == "0.4", ".40",simulation1_all_sum$tau)

#Plotting
simulation1_plot1 <- simulation1_all_sum %>% 
  filter(tau == ".25") %>% 
  plot_standard_by_name(
  title = "Baseline against heterogeneity of thresholds, tau = .25")

simulation1_plot2 <- simulation1_all_sum %>% 
  filter(name == "Heterogeneity of thresholds") %>% 
  plot_standard(title = "Heterogeneity of thresholds with varying tau")

randomtau_30 <- read.csv("randomtau_neu_rep100_tau30.csv")
randomtau_30 %>% filter(network > 25 & network < 37) %>% 
  ggplot(aes(round, adopters))+
  geom_line(size = 1.2)+
  theme_minimal()+
  theme(text = element_text(size = 15, family = "serif"), legend.key.height = unit(1,"cm"))+
  scale_color_brewer(palette = "PuOr")+
  labs(title = "Random tau 30", x = "Rounds", y = "Number of activated nodes")+
  facet_wrap(.~network)


plots_sim1 <- ggarrange(simulation1_plot1,simulation1_plot2, ncol = 2)
#

####Simulation 2 - stochastic thresholds####
#Loading file
simulation2_all_sum <- read.csv("simulation2_all_sum.csv")
simulation2_all_sum$tau <- ifelse(simulation2_all_sum$tau == "0.25", ".25",simulation2_all_sum$tau)
simulation2_all_sum$tau <- ifelse(simulation2_all_sum$tau == "0.3", ".30",simulation2_all_sum$tau)
simulation2_all_sum$tau <- ifelse(simulation2_all_sum$tau == "0.35", ".35",simulation2_all_sum$tau)
simulation2_all_sum$tau <- ifelse(simulation2_all_sum$tau == "0.4", ".40",simulation2_all_sum$tau)

#Plotting
simulation2_plot1 <- simulation2_all_sum %>% 
  filter(tau == ".25") %>% 
  plot_standard_by_name(title = "Simulation 2 - baseline against stochastic thresholds, tau = .25")+
  labs(title = str_wrap("Simulation 2 - baseline against stochastic thresholds, tau = .25",35))

simulation2_plot2 <- simulation2_all_sum %>% 
  filter(name == "Stochastic thresholds") %>% 
  plot_standard(title = "Simulation 2 - stochastic thresholds with varying tau")+
  labs(title = str_wrap("Simulation 2 - stochastic thresholds with varying tau",35))

plots_sim2 <- ggarrange(simulation2_plot1,simulation2_plot2, ncol = 2)
#


####Simulation 3 - heterogeneity of influence####
#Loading file
simulation3_all_sum <- read.csv("simulation3_all_sum.csv")
simulation3_all_sum$tau <- ifelse(simulation3_all_sum$tau == "0.25", ".25",simulation3_all_sum$tau)
simulation3_all_sum$tau <- ifelse(simulation3_all_sum$tau == "0.3", ".30",simulation3_all_sum$tau)
simulation3_all_sum$tau <- ifelse(simulation3_all_sum$tau == "0.35", ".35",simulation3_all_sum$tau)
simulation3_all_sum$tau <- ifelse(simulation3_all_sum$tau == "0.4", ".40",simulation3_all_sum$tau)

#Plotting
plot_sim3 <- simulation3_all_sum %>% 
  plot_standard_by_name(title = "Simulation 3 - baseline against heterogeneous influence, tau = .25")

#

####Simulation 4 - heterogeneity of degree####
#Loading file
#Baseline
simulation0_tau25 <- read.csv("basetau_neu_rep100_tau25.csv")
simulation0_tau25_sum <- sum_data(simulation0_tau25, name = "Baseline", tau = ".25")

simulation4_tau25 <- read.csv("basetau_scalefree_tau25_new.csv")
simulation4_tau25_sum <- sum_data(simulation4_tau25, name = "Heterogeneity of degree", tau = ".25")

simulation4_all_sum <- rbind(simulation0_tau25_sum, simulation4_tau25_sum)

#Plotting
plot_sim41 <- simulation4_all_sum %>% 
  plot_standard_by_name(title = "Simulation 4 - baseline against heterogeneous degree, tau = .25")+
  labs(title = str_wrap("Simulation 4 - baseline against heterogeneous degree, tau = .25",35))

plot_sim42 <- simulation4_all_sum %>% 
  filter(name == "Heterogeneity of degree") %>% 
  plot_standard(title = "Simulation 4 - heterogeneity of degree")+
  labs(title = str_wrap("Simulation 4 - heterogeneity of degree                          ",35))

plots_sim4 <- ggarrange(plot_sim41, plot_sim42, ncol = 2)


#Plot for showing differences in repetitions 

network11 <- simulation4_tau25 %>% filter(network == 11) %>% 
  ggplot(aes(round, adopters, color = "tau_type"))+
  geom_line(size = 1.2)+
  theme_minimal()+
  theme(text = element_text(size = 15, family = "serif"), legend.key.height = unit(1,"cm"))+
  scale_color_brewer(palette = "PuOr", name = "tau", labels = ".25")+
  labs(title = "Simulation 4, network 11", x = "Rounds", y = "Number of activated nodes")

network57 <- simulation4_tau25 %>% filter(network == 57) %>% 
  ggplot(aes(round, adopters, color = "tau_type"))+
  geom_line(size = 1.2)+
  theme_minimal()+
  theme(text = element_text(size = 15, family = "serif"), legend.key.height = unit(1,"cm"))+
  scale_color_brewer(palette = "PuOr", name = "tau", labels = ".25")+
  labs(title = "Simulation 4, network 57", x = "Rounds", y = "Number of activated nodes")

appendixx <- ggarrange(network11, network57)

#

####Simulation 5 - Interaction between stochastic thresholds and heterogeneity of thresholds####
#Loading file
simulation5_all_sum <- read.csv("simulation5_all_sum.csv")
simulation5_all_sum$tau <- ifelse(simulation5_all_sum$tau == "0.25", ".25",simulation5_all_sum$tau)
simulation5_all_sum$tau <- ifelse(simulation5_all_sum$tau == "0.3", ".30",simulation5_all_sum$tau)
simulation5_all_sum$tau <- ifelse(simulation5_all_sum$tau == "0.35", ".35",simulation5_all_sum$tau)
simulation5_all_sum$tau <- ifelse(simulation5_all_sum$tau == "0.4", ".40",simulation5_all_sum$tau)

#Plotting
simulation5_plot1 <- simulation5_all_sum %>% 
  filter(tau == ".25") %>% 
  plot_standard_by_name(title = "Simulation 5 - baseline against stochastic & heterogeneous thresholds, tau = .25")+
  labs(title = str_wrap("Simulation 5 - baseline against stochastic & heterogeneous thresholds, tau = .25", 45))

simulation5_plot2 <- simulation5_all_sum %>% 
  filter(name == "Stochastic & heterogeneous thresholds") %>% 
  plot_standard(title = "Simulation 5 - stochastic & heterogeneous thresholds with varying tau")+
  labs(title = str_wrap("Simulation 5 - stochastic & heterogeneous thresholds with varying tau",45))

plots_sim5 <- ggarrange(simulation5_plot1,simulation5_plot2, ncol = 2)

#

####Simulation 6 - Interaction between heterogeneity of thresholds and of influence####
#Loading file
simulation6_all_sum <- read.csv("simulation6_all_sum.csv")


#Plotting
plot_sim6 <- simulation6_all_sum %>% 
  plot_standard_by_name(title = "Simulation 6 - baseline against heterogeneity of degree and stochastic thresholds")

#
####Simulation 7 - Interaction between heterogeneity degree and heterogeneity of thresholds####
#Loading file
simulation7_all_sum <- read.csv("simulation7_all_sum.csv")

#Plotting
plot_sim7 <- simulation7_all_sum %>% 
  plot_standard_by_name(title = "Simulation 7 - baseline against heterogeneity of degree and of thresholds, tau = .25")

#


####All together####
#Loading all files 
#All simulations in one 

all_networks <- rbind(
  simulation0_tau25_sum, simulation0_tau30_sum, simulation0_tau35_sum, simulation0_tau40_sum, 
  simulation1_tau25_sum, simulation1_tau30_sum, simulation1_tau35_sum, simulation1_tau40_sum,
  simulation2_tau25_sum, simulation2_tau30_sum, simulation2_tau35_sum, simulation2_tau40_sum,
  simulation3_tau25_sum,
  simulation4_tau25_sum,
  simulation5_tau25_sum, simulation5_tau30_sum, simulation5_tau35_sum, simulation5_tau40_sum,
  simulation6_tau25_sum,
  simulation7_tau25_sum,
  simulation8_tau25_sum)

write.csv(all_networks, "all_data.csv")

#


####Old things####
all_simulations <- read.csv("all_networks.csv")
ggplot(all_simulations, aes(round, sumadopt, color = name))+
  geom_line(size = 1.2)+
  #geom_errorbar(stat = "summary", fun.data = mean_se)+
  geom_errorbar(aes(ymin=sumadopt-sd, ymax=sumadopt+sd), width=.2,
                position=position_dodge(.9), alpha = 0.25) +
  scale_colour_pander()+
  theme_minimal()+
  theme(legend.title = element_blank())+
  geom_hline(yintercept = 22500*0.5, linetype = 2, alpha = 0.7)+
  geom_text(aes(label ="50 % activation", x = 98, y = 12000), color = "Black", alpha = 0.01, size = 3)+
  geom_hline(yintercept = 22500*0.25, linetype = 2, alpha = 0.7)+
   geom_text(aes(label ="25 % activation", x = 98, y = 6300), color = "Black", alpha = 0.01, size = 3) +
  geom_hline(yintercept = 22500*0.75, linetype = 2, alpha = 0.7)+
   geom_text(aes(label ="75 % activation", x = 98, y = 17500),color = "Black", alpha = 0.01, size = 3) 

