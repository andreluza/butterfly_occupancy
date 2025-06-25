# -----------------------------------------------------------------------------

# Simulation study 1 - 1: replication of D&S Simulations higher smooth autocorrelation

# ----------------------------------------------------------------

# -------------------------------------------------------

# Bind the output of simulations Doser & Stoudt 2023
# (there were three breaks  during the processing, so now I need to bind the results).
#                   
# -------------------------------------------------------


rm(list = ls())
require(here)
require(abind)

# Load the data sets
load(file = here ("model_output", "output_simulations", "smooth_sims_D&S", "sim-data-correct.rda"))

# bind data
# Create new environments (three as there were 3 stops in the process)
env1 <- new.env()
env2 <- new.env()
#env3 <- new.env()

# Load the .rda files in each environment
load(file = here ("model_output", "output_simulations", "smooth_sims_D&S",
                  "sim-mixed-stPGOcc-results-SimSce1187.rda"), envir = env1)
load(file = here ("model_output", "output_simulations", "smooth_sims_D&S",
                  "sim-mixed-stPGOcc-results-SimSce1600.rda"), envir = env2)
#load(file = here ("model_output", "output_simulations", "sims_D&S",
#                  "sim-mixed-stPGOcc-results-SimSce1600.rda"), envir = env3)

# List objects in each environment
objects_env1 <- ls(env1)
objects_env1 <- objects_env1[-which(objects_env1 %in% c("alpha", "beta","scenario.vals"))]
objects_env2 <- ls(env2)
objects_env2 <- objects_env2[-which(objects_env2 %in% c("alpha", "beta","scenario.vals"))]
#objects_env3 <- ls(env3)
#objects_env3 <- objects_env3[-which(objects_env3 %in% c("alpha", "beta","scenario.vals"))]

# set the correct order of the scenarios (in the simulations I put lower phi first phi.vals <- (c(3 / 6, 3 / 3) ), whereas D&S did the converse phi.vals <- c(3 / .2, 3 / .8))

correct_order <- c(3,4,
                  1,2,
                  7,8,
                  5,6,
                  11,12,
                  9,10,
                  15,16,
                  13,14
                  )

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
    
    # get each object of environment 3
    #data3 <- get(objects_env3[i], envir = env3)
    # filter simulation that ended
    
    if (choose_dim == 1) {
    
        # filter simulation that ended
        data1 <- data1 [which(is.na(apply (data1,choose_dim,sum))!=T),,correct_order]
        data2 <- data2 [which(is.na(apply (data2,choose_dim,sum))!=T),,correct_order]
        #data3 <- data3 [which(is.na(apply (data3,choose_dim,sum))!=T),,]
        
        # Merged data
        merged_data <- abind(data1, data2, 
                             #data3, 
                             along = choose_dim )
        
    } else if (choose_dim == 2) {
      
      # filter simulation that ended
      data1 <- data1 [,which(is.na(apply (data1,choose_dim,sum))!=T),correct_order]
      data2 <- data2 [,which(is.na(apply (data2,choose_dim,sum))!=T),correct_order]
      #data3 <- data3 [,which(is.na(apply (data3,choose_dim,sum))!=T),]
      
      # Merged data
      merged_data <- abind(data1, data2,#data3, 
                           along = choose_dim )
      
    } else if (choose_dim == 3) {
      
      # filter simulation that ended
      data1 <- data1 [,,which(is.na(apply (data1,choose_dim,sum))!=T),correct_order]
      data2 <- data2 [,,which(is.na(apply (data2,choose_dim,sum))!=T),correct_order]
      #data3 <- data3 [,,which(is.na(apply (data3,choose_dim,sum))!=T),]
      
      # Merged data
      merged_data <- abind(data1, data2, #data3, 
                           along = choose_dim )
      
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
#rm(env3)
rm(objects_env1)
rm(objects_env2)
#rm(objects_env3)

ls()

# save 
save.image(file = here ("model_output", "output_simulations", "smooth_sims_D&S",
                                  "sim-mixed-stPGOcc-results-merged.rda"))

rm(list=ls())

