library(pacman)
p_load(tidyverse, network, igraph, intergraph, tidygraph, ggraph, ggplot2, ggthemes, comprehenr)

# Generate CA network with von Neumann neighborhoods 
generate_neumann <- function(n, p, nei){
  neumann <-  make_lattice(c(sqrt(n),sqrt(n)), nei = nei) %>% 
    rewire(each_edge(p = p)) 
  return(neumann)
}

# Create high connected nodes based on existing adjacency matrix 
create_connected_nodes <- function(adj, n_nodes, connectedness){
  ids = sample(nrow(adj), n_nodes) # Sample nodes to be highly connected
  adj_as_net = asNetwork(graph_from_adjacency_matrix(adj)) # make adj matrix into network 
  
  for (id in ids){
    n_new_connections = connectedness - length(get.neighborhood(adj_as_net, id)) # n_new_con = total_con - current_con
    new_connections = sample(nrow(adj), n_new_connections) # Create new connections for the node
    
    for (connection in new_connections){
      # Make the new connection in the matrix
      adj[id,connection] <- 1
      adj[connection,id] <- 1
      # Sample node to remove tie from 
      ones <- which(adj == 1, arr.ind = T)
      tie = sample(ones, 1)
      while (tie %in% ids){ # Make sure that tie to be removed is not high connected node
        tie = sample(ones, 1)
      }
      row = ones[tie,][[1]]
      col = ones[tie,][[2]]
      # Remove tie between two nodes
      adj[row,col] <- 0
      adj[col,row] <- 0
    }
  }
  return (list(adj, ids))
}

contagion_sim <- function(tau_type = "random_tau", # Takes values "random_tau" or "base_tau" 
                          n_seeds = 1, # Number of initial adopters
                          n_high = 0, #Number of high influence nodes, either status or degree
                         
                          high_status = F, # High status node with high influence
                          high_tau_perc = 1, # How influential high status nodes should be 
                          
                          high_degree = F, # High status nodes with high degree
                          connectedness = 0, # How connected high degree nodes should be
                          
                          rep = 50, # Number of repetitions - new network is generated every time
                          rounds = 50, # n rounds of the simulation 
                          tau = 0.33, # Threshold (mean of gaussian distribution or just threshold for all nodes)
                          n = 3600, # Number of nodes in the network - preferably a number with a natural square root
                          nei = 2, # degree of neighborhood connectedness --> 12 neighbors in von neumann 
                          p = 0.1){ # probability of rewiring
    
  # Ensure existence of data folder
  dir.create("data", showWarnings = F)
  
  # Prepare dataframe for degree distributions 
  degree_distribution <- data.frame()
  
  for (i in 1:rep){
    #generate network
    network <- generate_neumann(n, p, nei)
    networknetwork <- asNetwork(network)
   
    #### preparing the matrix ####
    adj <- as.matrix(as_adjacency_matrix(network, type = "both", sparse = T))
    diag(adj) <- 0 
    
    #### HGIH INFLUENCE NODES WITH HIGH DEGREE ####
    if (high_degree == T){
      mat_ids = create_connected_nodes(adj, n_high, connectedness)
      adj = mat_ids[1][[1]]
      networknetwork <- asNetwork(graph_from_adjacency_matrix(adj))
      id_connected_nodes <- mat_ids[2][[1]]
      
      #Preparing folder
      folder = file.path("data", paste("nHigh", n_high,
                                       "highDegree", high_degree, 
                                       "NConnectionsHigh", connectedness, 
                                       sep = "_"))
    }
    
    # Prepare neighbors and the percentage of the network that is activated 
    nei_activated_perc <- adj / rowSums(adj)
    
    #### High status nodes soccult way####
    if (high_status == T){
      # Influence of high status nodes can be a varying portion of mean tau 
      tau_high <- tau * high_tau_perc
      extra_influ <- ((tau_high * 12)-1)*n_high
      penalty <- extra_influ / n /12
      nei_activated_perc <- ifelse(nei_activated_perc !=0, nei_activated_perc - penalty, nei_activated_perc) 
      #now everybody has lost the amount of influence that is distrubted to high status nodes
      
      #defining highstat nodes
      high_status_nodes <- sample(seq_len(n), size = n_high)
      #looping through all high stat nodes and assigning new value to them
      for (b in high_status_nodes){
        node <- b
        neighbors <- get.neighborhood(networknetwork, node)
        nei_activated_perc[neighbors, node] <- tau_high #now all high nodes have tau as influence
      }
      
      ## Preparing folder
      folder = file.path("data", paste("nHigh", n_high,
                                       "highStatus", high_status, 
                                       "HighNodeTauPerc", high_tau_perc,
                                       sep = "_"))
    }
    ## Creating folder for baseline simulations
    if (high_status == F & high_degree == F){
      folder = file.path("data", paste("NoHigh", "tau", tau, "nSeeds", n_seeds, sep = "_"))
    }
    # Creating folder for storing data 
    dir.create(folder, showWarnings = F)
    
    # Saving degree distribution for this iteration 
    degree = data.frame(node=rep(0,n), neighbors=rep(0,n), n_neighbors = rep(0,n))
    
    for (j in 1:n){
      neighbors = get.neighborhood(networknetwork, j)
      degree$node[j] = j
      degree$neighbors[j] = toString(neighbors)
      degree$n_neighbors[j] = length(neighbors)
    }
    distribution <- degree %>% 
      dplyr::group_by(n_neighbors) %>% 
      dplyr::summarise(n=n())
    distribution$repetition <- i
    degree_distribution <- rbind(degree_distribution, distribution)
    
    #### INITIAL ADOPTER #### 
    adopters <- rep(F, n)
    # Choose n_seeds initial adopters at random
    initial_adopters <- base::sample(seq_len(n), size = n_seeds) 
  
    for (adopter in initial_adopters){
      initial_neighbors <- get.neighborhood(networknetwork, adopter) #get neighbors for each initial adopter
      adopters[c(adopter, initial_neighbors)] <- T # Set neighbors as adopters
    }
    print(paste("Starting iteration number", i, "at time", Sys.time()))
    
    
    #### PREPARE LIST FOR STORING ADOPTERS EACH ROUND ####
    adopt <- vector(mode = "list", length = rounds) 
    adopt[[1]] <- adopters
    
    #### Create tau vector 
    
    if(tau_type == "base_tau"){
      #preparing tau list
      tau_vec <- as.vector(rep(tau, nrow(adj)))
    }
    if(tau_type == "random_tau"){ 
      tau_norm <- rnorm(nrow(adj), tau, 0.16)
      tau_vec<- as.vector(sample(tau_norm))
    }
    
    #### Run simulation for rounds rounds
    for (t in 2:rounds) {
      adopt[[t]] <- ifelse(adopters, TRUE, ((nei_activated_perc %*% adopt[[t - 1]]) >= tau_vec)) 
      
      # Save information about iteration i, round t
      df <- data.frame(
        network = i,
        round = t,
        adopters = sum(adopt[[t]]),
        high_status = high_status,
        high_degree = high_degree,
        tau_type = paste(tau_type)
      )
      
      if (t == 2) {
        adopt_dat <- df
      }else{
        adopt_dat <- rbind(adopt_dat, df)
      }
    }
    
    # Save information about iteration i, round 1
    df_1 <- data.frame(
      network = i,
      round = 1,
      adopters = sum(adopt[[1]]),
      high_status = high_status,
      high_degree = high_degree,
      tau_type = paste(tau_type)
    )
    adopt_dat <- rbind(df_1, adopt_dat)
    if (i == 1){
      all_adopters <- adopt_dat
    }else{
      all_adopters <- rbind(all_adopters, adopt_dat)
    }
    
  }
  # Save simulation results (diffusion) and degree distributions
  write.csv(all_adopters, file.path(folder, "simulation_results.csv"))
  write.csv(degree_distribution, file.path(folder, "degree_distribution.csv"))
  return(all_adopters)
}
