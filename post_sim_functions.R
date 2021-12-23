# Functions for summarizing dataframes, calculating point estimates and plotting

# This function is for plotting a single simulation
plot_standard <- function(summed_dataframe, title = "",
                          x_name = "Rounds", y_name = "Number of activated nodes"){
  data = summed_dataframe
  plot <- ggplot(data, aes(round, sumadopt, color = str_wrap(tau,20)))+
    geom_line(size = 1.2)+
    theme_minimal()+
    theme(text = element_text(size = 15, family = "serif"), legend.key.height = unit(1,"cm"))+
    scale_color_brewer(palette = "PuOr")+
    labs(title = title, x = x_name, y = y_name) # , color = "tau"
  return(plot)
}

# This function is for plotting several simulation by the name given in summing. 
plot_standard_by_name <- function(summed_dataframe, title = "",
                                  x_name = "Rounds", y_name = "Number of activated nodes"){
  data = summed_dataframe
  plot <- ggplot(data, aes(round, sumadopt, color = str_wrap(name,20)))+
    geom_line(size = 1.2)+
    theme_minimal()+
    theme(text = element_text(size = 15, family = "serif"), legend.key.height = unit(1,"cm"))+
    scale_color_brewer(palette = "Dark2")+
    labs(title = title, x = x_name, y = y_name, color = "Simulation")
  return(plot)
}





# Plotting degree distribution
plot_degree_distribution <- function(mean_df, title, xlab="Number of neighbors", ylab="Count"){
  plot <- ggplot(mean_df) +
            geom_bar(aes(x= n_neighbors,y=mean_n), stat="identity", fill ="white", color="black", width=1) +
            geom_errorbar(aes(x=n_neighbors, ymin=mean_n-sd_n, ymax=mean_n+sd_n), width=0.7, color="gray20")+
            labs(x=xlab, y = ylab, title=title)+
            theme_calc()+
            theme(text = element_text(size = 15, family = "serif"))+
            scale_x_continuous(breaks = seq(0,length(mean_df$n_neighbors)), minor_breaks = 0)
  return(plot)
  }


# This function is for summing the data, name is what it will be plot by 
sum_data <- function(dataframe, name = "", baseline, status, degree){
  summed_data <- dataframe %>% group_by(round) %>% summarise(sumadopt = mean(adopters),
                                                             sd = sd(adopters),
                                                             name = name,
                                                             baseline = baseline, 
                                                             status = status, 
                                                             degree = degree)
  return(summed_data)
}

# Automatically sum simulation results, takes list of folders as input
# Outputs the summed results and the degree distributions 
sum_multiple_results <- function(folders){
  all_data_sum = data.frame()
  degree_distributions <- data.frame()
  for (folder in folders){
    # Split folder name to extract information 
    split = strsplit(folder, "_|/") # Split by _ OR /
    results = read.csv(paste(folder, "/simulation_results.csv", sep=""))
    degree_files = list.files(folder, glob2rx("degree_*.csv"))
    
    degree_distribution <- data.frame()
    for (file in degree_files){
      degree = read.csv(paste(folder, file,sep="/"))
      # print(head(degree))
      degree <- degree[,2:3]
      degree_distribution <- rbind(degree_distribution, degree)
    }
    degree_distribution$baseline <- FALSE
    degree_distribution$status <- FALSE
    degree_distribution$degree <- FALSE
    degree_distribution$name <- paste(split[[1]][length(split[[1]])], "seeds")
    # Baseline files
    if (split[[1]][2] == "NoHigh"){
      summed_data <- sum_data(results, name = paste(split[[1]][length(split[[1]])], "seeds"), T, F, F)
      all_data_sum <- rbind(all_data_sum, summed_data)
      degree_distribution$baseline <- TRUE
      degree_distributions <- rbind(degree_distributions, degree_distribution)
    } #Finish summing baseline
    
    # High status files
    else if (split[[1]][4] == "highStatus"){
      summed_data <- sum_data(results, name = paste(split[[1]][3], 
                                                    "seeds"), 
                              F, split[[1]][length(split[[1]])], F)
      all_data_sum <- rbind(all_data_sum, summed_data)
      degree_distribution$status <- split[[1]][length(split[[1]])]
      degree_distributions <- rbind(degree_distributions, degree_distribution)
    } #Finish summing high status
    
    # High degree files 
    else if (split[[1]][4] == "highDegree"){
      summed_data <- sum_data(results, name = paste(split[[1]][3], 
                                                    "seeds"), 
                              F, F, split[[1]][length(split[[1]])])
      all_data_sum <- rbind(all_data_sum, summed_data)
      degree_distribution$degree <- split[[1]][length(split[[1]])]
      degree_distributions <- rbind(degree_distributions, degree_distribution)
    } #Finish summing high degree 
    
  } # Looped through all folder
  return (list(all_data_sum, degree_distributions)) #
}



