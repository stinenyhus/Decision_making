# Decision Making

## Project organization 
This repository contains code for running, analyzing and plotting simulation for our decision making exam. 
It is structured as follows

```
├── .gitignore                      # A list of files not uploaded to git
├── README.md                       # The top-level README for this project.
├── data                            # Folder for storing simulation outputs - automatically created when running a simulation
|   ├── NoHigh_tau_0.33_nSeeds_2    # Folder containing data from baseline simulation with 2 initial seeds 
|       ├── degree_distribution.csv # Degree distributions from all iterations
|       └── simulation_results.csv  # Propagation information for all iterations 
|   └── ...                         # total of n_simulations folders, automatically created
├── functions                       # The main folder for scripts
|   └── simulation.r                # Support functions and main functions for running simulations
|   └── post_sim_functions.R        # Functions for summarizing and plotting data after simulations are run 
├── run_simulations.Rmd             # Markdown file for running the simulations used in this study
└── sum_and_plot.Rmd                # Markdown file for creating plots and summary statistics of the simulations in this study
```

## Running the simulations
To run the simulations, run function contagion_sim with desired variables. 
The variables of the function contagion_sim that are held constant throughout all simulations in this study are 
- tau_type = "random_tau" (nodes are randomly assigned a threshold drawn from a normal distribution)
- rep = 50 (all conditions are run in 50 iterations that are then averaged)
- rounds = 50 (all simulations run for 50 rounds no matter how slowly/quickly the iteration runs)
- tau = 0.33 (since tau_type = "random_tau", this value refers to the mean of normal distribution. Standard deviation is hardcoded = 0.16)
- n = 3600 (network size is 3600 nodes for all simulations. This created 60 x 60 two-dimentional CA network)
- nei = 2 (most nodes will initially have 12 neighbors)
- p = 0.1 (rewiring probability is 10% for all simulations)
- n_seeds = 1 (the contagion is seeded with one node and its neighborhood)

### High degree 
High influence nodes with high degree are implemented by setting high_degree = T in contagion_sim and also providing a non-zero values for n_high (indicating number of high degree nodes in the network) and connectedness (indicating how many connections high degree nodes should have). 

### High status
High influence nodes with high status are implemented by setting high_status = T in contagion_sim and also providing a non-zero values for n_high (indicating number of high status nodes in the network) and high_tau_perc (indicating how much influence out of mean tau high status should have. high_tau_perc = 1 means that high status nodes have 100% of mean tau as influence). 

### Post hoc
Additionally to the main simulations, post hoc simulations are run. These hold influence and degree constant at 100% of tau and 24 connections, respectively, and increase number of high status nodes to 2, 4 and 8 %. 