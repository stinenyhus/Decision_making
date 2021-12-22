library(pacman)
p_load(comprehenr)

generate_neumann <- function(n, p, nei){
  neumann <-  make_lattice(c(sqrt(n),sqrt(n)), nei = nei) %>% 
    rewire(each_edge(p = p)) 
  return(neumann)
}

generate_scalefree <- function(n, degrees, gamma, p = 0){
  degs <- sample(1:degrees, n, replace=TRUE, prob=(1:degrees)^-gamma)
  if (sum(degs) %% 2 != 0) { degs[1] <- degs[1] + 1 } ### Note, that we correct the degree sequence if its sum is odd
 # degs <- degs * 2
  scale_free <- sample_degseq(degs, method = "vl") 
  scale_free <- scale_free %>% rewire(each_edge(p = p)) 
  return(scale_free)
} 

create_connected_nodes <- function(adj, n_nodes, connectedness){
  ids = sample(nrow(adj), n_nodes) #Sample nodes to be highly connected
  adj_as_net = asNetwork(graph_from_adjacency_matrix(adj)) #make adj matrix into network 
  for (id in ids){
    n_new_connections = connectedness - length(get.neighborhood(adj_as_net, id)) # n_new_con = total_con - current_con
    new_connections = sample(nrow(adj), n_new_connections) #Create new connections for the node
    print(paste("Number of sampled new connections is" , length(new_connections)))
    #Check whether the new connections are already neighbors - DELETE THIS PART? SEEMS UNNECESSARY 
    connections = comprehenr::to_vec(for(n in new_connections) if (adj[id,n]==1) n+1 else n)
    print(paste("Number of non-neighbors is new connections is", length(connections)))
    while (new_connections != connections){ #Update until all connections are new
      connections = comprehenr::to_vec(for(n in connections) if (adj[id,n]==1) n+1 else n)
      print(paste("Number of non-neighbors is new connections is", length(connections)))
    }
    for (connection in connections){
      # Make the new connection in the matrix
      adj[id,connection] <- 1
      adj[connection,id] <- 1
      # Sample node to remove tie from 
      ones <- which(adj == 1, arr.ind = T)
      tie = sample(ones, 1)
      while (tie %in% ids){ #Make sure that tie to be removed is not high connected node
        tie = sample(ones, 1)
      }
      row = ones[tie,][[1]]
      col = ones[tie,][[2]]
      # Remove tie between two nodes
      adj[row,col] <- 0
      adj[col,row] <- 0
    }
  }
  return (adj)
}

contagion_sim <- function(tau_type = "base_tau", # Takes values "random_tau" or "base_tau" 
                          high_node = F, # High status node old way 
                          n_connected_nodes = 0, #Number of highly connected nodes
                          connectedness = 0, 
                          rep, # Number of repetitions - new network is generated every time
                          net_type = "neumann", # Takes value "neumann" or "scale_free"
                          rounds, # n rounds of the simulation 
                          tau = 0, # Threshold (mean of gaussian distribution or just threshold for all nodes)
                          n, # Number of nodes in the network - preferably a number with a natural square root
                          nei, # degree of neighborhood connectedness
                          p = 0, # probability of rewiring - keep zero for new high status node?
                          degrees = 0, # Only use for scale free network 
                          gamma = 0, # Only use for scale free network
                          stochastic = F){ # whether or not thresholds should be stochastic 
  for (i in 1:rep){
    #generate network
    if(net_type == "neumann"){
      network <- generate_neumann(n, p, nei)
      networknetwork <- asNetwork(network)
    }
    if(net_type == "scale_free"){
      network <- generate_scalefree(n, degrees, gamma, p)
      networknetwork <- asNetwork(network)
    }
    
    #### preparing the matrix ####
    adj <- as.matrix(as_adjacency_matrix(network, type = "both", sparse = T))
    diag(adj) <- 0 
    
    #### HIGHLY CONNECTED NODES DECISION MAKING WAY ####
    if (n_connected_nodes != 0){
      adj = create_connected_nodes(adj, n_connected_nodes, connectedness)
    }
    nei_activated_perc <- adj / rowSums(adj)
    
    
    #### INITIAL ADOPTER #### 
    adopters <- rep(F, n)
    # Choose a person at random
    initial_adopter <- base::sample(seq_len(n), size = 1) #takes one node from the sequence of nodes and classifies it as the inital 
    #adopter
    initial_neighbors <- get.neighborhood(networknetwork, initial_adopter) #shows the 3 neighbors that the infected adopter has
    adopters[c(initial_adopter, initial_neighbors)] <- T #setting all neighbors to the initial infected node as infected
    print(paste("Initial adopters found for repetition number ", i))
    
    
    #### High status nodes soccult way####
    if (high_node == T){
      n_high <- n/12
      extra_influ <- ((tau * 12)-1)*n_high
      penalty <- extra_influ / n /12
      nei_activated_perc <- ifelse(nei_activated_perc !=0, nei_activated_perc - penalty, nei_activated_perc) 
      #now everybody has lost the amount of influence that is distrubted to high status nodes now
      #defining highstat nodes
      high_status_nodes <- sample(seq_len(n), size = n_high)
      #looping through all high stat nodes and assigning new value to them
      for (b in high_status_nodes){
        node <- b
        neighbors <- get.neighborhood(networknetwork, node)
        nei_activated_perc[neighbors, node] <- tau #now all high nodes have tau as influence
      }
    }
    
    
    #### PREPARE LIST FOR STORING ADOPTERS EACH ROUND ####
    adopt <- vector(mode = "list", length = rounds) 
    adopt[[1]] <- adopters
    
    ###### SIMULATIONS UNDER DIFFERENT CONDITIONS ######
    
    #### BASELINE NETWORK ####
    if(tau_type == "base_tau" & stochastic == F){
      #preparing tau list
      tau_vec <- as.vector(rep(tau, nrow(adj)))
      for (t in 2:rounds) {
        print(t)
        adopt[[t]] <- ifelse(adopters, TRUE, ((nei_activated_perc %*% adopt[[t - 1]]) >= tau_vec)) 
        
        df <- data.frame(
          network = i,
          round = t,
          adopters = sum(adopt[[t]]),
          high_node = high_node,
          stochastic = paste(stochastic),
          tau_type = paste(tau_type)
        )
        
        if (t == 2) {
          adopt_dat <- df
        }else{
          adopt_dat <- rbind(adopt_dat, df)
        }
      }
      df_1 <- data.frame(
        network = i,
        round = 1,
        adopters = sum(adopt[[1]]),
        high_node = high_node,
        stochastic = paste(stochastic),
        tau_type = paste(tau_type)
      )
      adopt_dat <- rbind(df_1, adopt_dat)
      if (i == 1){
        all_adopters <- adopt_dat
      }else{
        all_adopters <- rbind(all_adopters, adopt_dat)
      }
    }
    
    #### NORMALLY DISTRIBUTED THRESHOLDS ####
    if(tau_type == "random_tau" & stochastic == F){ 
      tau_norm <- rnorm(nrow(adj),tau, tau/3)
      tau_vec<- as.vector(sample(tau_norm))
      for (t in 2:rounds) {
        adopt[[t]] <- ifelse(adopters, TRUE, ((nei_activated_perc %*% adopt[[t - 1]]) >= tau_vec)) 
        
        df <- data.frame(
          network = i,
          round = t,
          adopters = sum(adopt[[t]]),
          high_node = high_node,
          stochastic = paste(stochastic),
          tau_type = paste(tau_type)
        )
        if (t == 2) {
          adopt_dat <- df
        }else{
          adopt_dat <- rbind(adopt_dat, df)
        }
      }
      df_1 <- data.frame(
        network = i,
        round = 1,
        adopters = sum(adopt[[1]]),
        high_node = high_node,
        stochastic = paste(stochastic),
        tau_type = paste(tau_type)
      )
      adopt_dat <- rbind(df_1, adopt_dat)
      if (i == 1){
        all_adopters <- adopt_dat
      }else{
        all_adopters <- rbind(all_adopters, adopt_dat)
      }
    }
    
    #### STOCHASTIC THRESHOLD - UNIFORM OR NORM DISTRIB ####
    if(stochastic == T){
      if (tau_type == "base_tau"){
        tau_vec = as.vector(rep(tau, nrow(adj)))
      }
      else if (tau_type == "random_tau"){
        tau_vec = rnorm(nrow(adj),tau, tau/3) #Some taus are negative which is not good! 
        #They should only be pos so therefore sd = tau/3
      }
      for (t in 2:rounds){
        print(t)
        list <- adopt[[t-1]]
        m = nei_activated_perc %*% adopt[[t - 1]]
        q = m - tau_vec + 0.5
        for (s in 1:length(list)){
          L = 1/(1+exp((.5-q[s])*10))
          activation = rbinom(1,1,prob = L)
          if (activation == 1) {
            list[s] = TRUE}
          else if (activation == 0) {#If we remove this else if statement, maybe they cant turn of?
            list[s] = FALSE}
        }
        adopt[[t]] <- list
        
        df <- data.frame(
          network = i,
          round = t,
          adopters = sum(adopt[[t]]),
          high_node = high_node,
          stochastic = paste(stochastic),
          tau_type = paste(tau_type)
        )
        if(t ==2){
          adopt_dat <- df
        }else{
          adopt_dat <- rbind(adopt_dat, df)
        }
      }
      df_1 <- data.frame(
        network = i,
        round = 1,
        adopters = sum(adopt[[1]]),
        high_node = high_node,
        stochastic = paste(stochastic),
        tau_type = paste(tau_type)
      )
      adopt_dat <- rbind(df_1, adopt_dat)
      
      if(i ==1){
        all_adopters <- adopt_dat
      } else {
        all_adopters <- rbind(all_adopters, adopt_dat)
      }
      
    }
    
  }
  return(all_adopters)
}



plot_standard <- function(summed_dataframe, title = "",
                          x_name = "Rounds", y_name = "Number of activated nodes"){
    data = summed_dataframe
    plot <- ggplot(data, aes(round, sumadopt, color = str_wrap(tau,20)))+
    geom_line(size = 1.2)+
    theme_minimal()+
    theme(text = element_text(size = 15, family = "serif"), legend.key.height = unit(1,"cm"))+
    scale_color_brewer(palette = "PuOr")+
    labs(title = title, x = x_name, y = y_name, color = "tau")
  return(plot)
}

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


sum_data <- function(dataframe, name = "", tau = ""){
  summed_data <- dataframe %>% group_by(round) %>% summarise(sumadopt = mean(adopters),
                                                             sd = sd(adopters),
                                                             name = name,
                                                             tau = tau)
  return(summed_data)
}


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


