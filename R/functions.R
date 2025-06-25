

# -----------------------------------------------
# function to find closest values (useful to consider minimum N observers (target) when setting NAs)
# keep consistency on observations
find_closest_values <- function(values, target, n = 50) {
  # Ensure values is a vector
  #if (!is.vector(values)) stop("The input 'values' must be a vector")
  
  # Calculate the absolute difference between each value and the target
  differences <- abs(values - target)
  
  # Order the values by the difference
  ordered_indices <- order(differences)
  
  # Select the first n indices
  closest_indices <- ordered_indices[1:n]
  
  # Return the corresponding values
  closest_values <- values[closest_indices]
  
  return(closest_indices)
}

# -----------------------------------------------------------
#--------------using sampling library rather than base------------------
# corrected sampling function: source: 
#https://www.r-bloggers.com/2024/08/sampling-without-replacement-with-unequal-probabilities-by-ellis2013nz/



library(sampling)
sample_brewer <- function(x, size, prob, replace = FALSE, keep = FALSE){
  pik <- prob / sum(prob) * size
  s <- UPbrewer(pik)
  the_sample <- x[which(s == 1)]
  return(the_sample)
}


# probability mass of a Poisson distribution
function_prob_poisson <- function (lambda,k) (lambda^k)*exp(-lambda)/factorial(k)


# gaussian function for phenology
seasonal_effect <- function(j, peak, spread) { 
  scaling_factor * exp(-decay * ((j - peak) / spread)^2) 
}

