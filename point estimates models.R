library(pacman)
p_load(tidyverse)
source("useful_function_soccult.r")

####Simulation 0 - baseline####
simulation0_25 <- read.csv("basetau_neu_rep100_tau25.csv")
simulation0_25_summed <- sum_data(simulation0_25, name = "Sim 0, tau.25")
simulation0_25_estimates <- calculate_point_estimates(simulation0_25, simulation0_25_summed)
simulation0_30 <- read.csv("basetau_neu_rep100_tau30.csv")
simulation0_30_summed <- sum_data(simulation0_30, name = "Sim 0, tau.30")
simulation0_30_estimates <- calculate_point_estimates(simulation0_30, simulation0_30_summed)
simulation0_35 <- read.csv("basetau_neu_rep100_tau35.csv")
simulation0_35_summed <- sum_data(simulation0_35, name = "Sim 0, tau.35")
simulation0_35_estimates <- calculate_point_estimates(simulation0_35, simulation0_35_summed)
simulation0_40 <- read.csv("basetau_neu_rep100_tau40.csv")
simulation0_40_summed <- sum_data(simulation0_40, name = "Sim 0, tau.40")
simulation0_40_estimates <- calculate_point_estimates(simulation0_40, simulation0_40_summed)
#

####Simulation 1 - heterogeneity of thresholds####
simulation1_25 <- read.csv("randomtau_neu_rep100_tau25.csv")
simulation1_25_summed <- sum_data(simulation1_25, name = "Sim 1 tau 25")
simulation1_25_estimates <- calculate_point_estimates(simulation1_25, simulation1_25_summed)
simulation1_30 <- read.csv("randomtau_neu_rep100_tau30.csv")
simulation1_30_summed <- sum_data(simulation1_30, name = "Sim 1 tau 30")
simulation1_30_estimates <- calculate_point_estimates(simulation1_30, simulation1_30_summed)
simulation1_35 <- read.csv("randomtau_neu_rep100_tau35.csv")
simulation1_35_summed <- sum_data(simulation1_35, name = "Sim 1 tau 35")
simulation1_35_estimates <- calculate_point_estimates(simulation1_35, simulation1_35_summed)
simulation1_40 <- read.csv("randomtau_neu_rep100_tau40.csv")
simulation1_40_summed <- sum_data(simulation1_40, name = "Sim 1 tau 40")
simulation1_40_estimates <- calculate_point_estimates(simulation1_40, simulation1_40_summed)
#
22500*0.25

####Simulation 2 - stochastic thresholds####
simulation2_25 <- read.csv("basetau_neu_rep100_stochastic_tau_25.csv")
simulation2_25_summed <- sum_data(simulation2_25, name = "Sim 2 tau 25")
simulation2_25_estimates <- calculate_point_estimates(simulation2_25, simulation2_25_summed)
simulation2_30 <- read.csv("basetau_neu_rep100_stochastic_tau30.csv")
simulation2_30_summed <- sum_data(simulation2_30, name = "Sim 2 tau 30")
simulation2_30_estimates <- calculate_point_estimates(simulation2_30, simulation2_30_summed)
simulation2_35 <- read.csv("basetau_neu_rep100_stochastic_tau_35.csv")
simulation2_35_summed <- sum_data(simulation2_35, name = "Sim 2 tau 35")
simulation2_35_estimates <- calculate_point_estimates(simulation2_35, simulation2_35_summed)
simulation2_40 <- read.csv("basetau_neu_rep100_stochastic_tau40.csv")
simulation2_40_summed <- sum_data(simulation2_25, name = "Sim 2 tau 40")
simulation2_40_estimates <- calculate_point_estimates(simulation2_40, simulation2_40_summed)
#

####Simulation 3 - heterogeneity of influence####
simulation3 <- read.csv("basetau_neu_rep100_highstatus_tau25.csv")
vec <- simulation3$adopters[simulation3$round == 100] > 100
trues <- which(vec, T)
simulation3_success <- simulation3[simulation3$network %in% trues,]
simulation3_success_summed <- sum_data(simulation3_success, name = "Sim 3")
simulation3_estimates_success <- calculate_point_estimates(simulation3_success, simulation3_success_summed)

#


####Simulation 4 - heterogeneity of degree####
simulation4 <- read.csv("basetau_scalefree_rep100_tau25.csv")
simulation4_summed <- sum_data(simulation4, name = "Sim 4")
simulation4_estimates <- calculate_point_estimates(simulation4, simulation4_summed)
#

####Simulation 5 - Interaction between stochastic thresholds and heterogeneity of thresholds####
simulation5_25 <- read.csv("randomtau_neu_rep100_stochastic_25.csv")
simulation5_25_summed <- sum_data(simulation5_25, name = "Sim 5 tau 25")
simulation5_25_estimates <- calculate_point_estimates(simulation5_25, simulation5_25_summed)
simulation5_30 <- read.csv("randomtau_neu_rep100_stochastic_tau30.csv")
simulation5_30_summed <- sum_data(simulation5_30, name = "Sim 5 tau 30")
simulation5_30_estimates <- calculate_point_estimates(simulation5_30, simulation5_30_summed)
simulation5_35 <- read.csv("randomtau_neu_rep100_stochastic_tau35.csv")
simulation5_35_summed <- sum_data(simulation5_35, name = "Sim 5 tau 35")
simulation5_35_estimates <- calculate_point_estimates(simulation5_35, simulation5_35_summed)
simulation5_40 <- read.csv("randomtau_neu_rep100_tau40_stochastic.csv")
simulation5_40_summed <- sum_data(simulation5_40, name = "Sim 5 tau 40")
simulation5_40_estimates <- calculate_point_estimates(simulation5_40, simulation5_40_summed)

#
####Simulation 6 - Interaction between heterogeneous thresholds and high status ####
simulation6 <- read.csv("randomtau_neu_rep100_highstatus.csv")
simulation6_summed <- sum_data(simulation6, name = "Sim 6")
simulation6_estimates <- calculate_point_estimates(simulation6, simulation6_summed)

#
####Simulation 7 - Interaction between heterogeneity degree and heterogeneity of thresholds####
simulation7 <- read.csv("randomtau_scalefree_rep100_tau25.csv")
simulation7_summed <- sum_data(simulation7, name = "Sim 7")
simulation7_estimates <- calculate_point_estimates(simulation7, simulation7_summed)

#
####Simulation 8 - Interaction between stochastic thresholds and heterogeneity of influence####
simulation8 <- read.csv("basetau_neu_rep100_highstatus_stochastic.csv")
simulation8_summed <- sum_data(simulation8, name = "Int. stochastic thresholds & heterogeneity of influence")
simulation8_estimates <- calculate_point_estimates(simulation8, simulation8_summed)

#
####All estimates####
all_estimates <- rbind(
  simulation0_25_estimates, simulation0_30_estimates, simulation0_35_estimates, simulation0_40_estimates,
  simulation1_25_estimates, simulation1_30_estimates, simulation1_35_estimates, simulation1_40_estimates,
  simulation2_25_estimates, simulation2_30_estimates, simulation2_35_estimates, simulation2_40_estimates,
  simulation3_estimates_success, 
  simulation4_estimates,
  simulation5_25_estimates, simulation5_30_estimates, simulation5_35_estimates, simulation5_40_estimates,
  simulation6_estimates,
  simulation7_estimates)
