# -----------------------------------------------------------------------------

# Simulation study 2: scenario 2

# -------------------------------------------------------------------------------
# -------------------------------------------------------

# Bind the output of simulations
# (there were three breaks  during the processing, so now I need to bind the results).
#                   
# -------------------------------------------------------

rm(list = ls())
require(here)
require(abind)

# Load the data sets
load(file = here ("model_output", 
                  "output_simulations", 
                  "sims_D&S", 
                  "sim-data-correct.rda"))

# bind data
# Create new environments (three as there were 3 stops in the process)
env1 <- new.env()
env2 <- new.env()
env3 <- new.env()

# Load the .rda files in each environment
load(file = here ("model_output", "output_simulations", "scenario_zero",
                  "sim-mixed-stPGOcc-results_1294.rda"), envir = env1)
load(file = here ("model_output", "output_simulations", "scenario_zero",
                  "sim-mixed-stPGOcc-results_1600.rda"), envir = env2)

# List objects in each environment
objects_env1 <- ls(env1)
objects_env1 <- objects_env1[-which(objects_env1 %in% c("alpha", "beta","scenario.vals"))]
objects_env2 <- ls(env2)
objects_env2 <- objects_env2[-which(objects_env2 %in% c("alpha", "beta","scenario.vals"))]

# Merge object
bind_DS_output <- lapply (seq (1,length(objects_env1)), function (i) {
    
  tryCatch ({
    
    # get each object of environment 1
    data1 <- get(objects_env1[i], envir = env1)
    # choose the dimension to filter (dimension of simulations)
    if (grepl("\\beta\\b",objects_env1[i])) {
      choose_dim <- 1
    } else {
      choose_dim <- max(which (dim(data1) == 100)) 
    }
    
    # get each object of environment 2
    data2 <- get(objects_env2[i], envir = env2)
    
    # filter simulation that ended
    
    if (choose_dim == 1) {
    
      # some simulations stopped in the intermediary scenarios, so we need to discard them and bind the complete results for all scenarios within one simulation
      # data1 [which(is.na(apply (data2,choose_dim,sum,na.rm=F))!=T)[1],,] # example
      data1 [which(is.na(apply (data2,choose_dim,sum,na.rm=F))!=T)[1],,] <- NA # set NA
      # data1[,,7] # check
      
      # filter simulation that ended
      data1 <- data1 [which(is.na(apply (data1,choose_dim,sum))!=T),,]
      data2 <- data2 [which(is.na(apply (data2,choose_dim,sum))!=T),,]
        
      # Merged data
      merged_data <- abind(data1, 
                           data2, 
                           along = choose_dim )
        
    } else if (choose_dim == 2) {
      
      # some simulations stopped in the intermediary scenarios, so we need to discard them and bind the complete results for all scenarios within one simulation
      # data1 [,which(is.na(apply (data2,choose_dim,sum,na.rm=F))!=T)[1],] # example
      data1 [,which(is.na(apply (data2,choose_dim,sum,na.rm=F))!=T)[1],] <- NA # set NA
      # data1[,7,] # check
      
      # choose data to bind
      data1 <- data1 [,which(is.na(apply (data1,choose_dim,sum))!=T),]
      data2 <- data2 [,which(is.na(apply (data2,choose_dim,sum))!=T),]
      
      # Merged data
      merged_data <- abind(data1, data2,
                           along = choose_dim )
      
    } else if (choose_dim == 3) {
      
      # some simulations stopped in the intermediary scenarios, so we need to discard them and bind the complete results for all scenarios within one simulation
      # data1 [,,which(is.na(apply (data2,choose_dim,sum,na.rm=F))!=T)[1],13] # example
      data1 [,,which(is.na(apply (data2,choose_dim,sum,na.rm=F))!=T)[1],] <- NA # set NA
      # data1[,,7,] # check
      
      # filter simulation that ended
      data1 <- data1 [,,which(is.na(apply (data1,choose_dim,sum))!=T),]
      data2 <- data2 [,,which(is.na(apply (data2,choose_dim,sum))!=T),]
      
      # Merged data
      merged_data <- abind(data1, data2, 
                           along = choose_dim )
      # check
      # tail(merged_data[,,7,])
      
    }
  
    merged_data
    
    
  }, error = function (e) return (NULL)) # check if any error comes out
    
})

# set names
names(bind_DS_output) <- objects_env1

# release the objects in the global environment
list2env(bind_DS_output,
         envir = .GlobalEnv)

# remove the big list to save memory
rm(bind_DS_output)
rm(dat)
rm(env1)
rm(env2)
rm(objects_env1)
rm(objects_env2)

# save 
save.image(file = here ("model_output", "output_simulations", "scenario_zero",
                                  "sim-mixed-stPGOcc-results-merged.rda"))

rm(list=ls())

