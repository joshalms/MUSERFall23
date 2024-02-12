
#Approach:
  #Increment viral load for all E and I (should be very short)



viral_load <- function(dat, at) {

  ## Attributes ##

  active <- get_attr(dat, "active") #list of all people who are doing something at this time step
  status <- get_attr(dat, "status") #List of all statuses with indices corresponding to IDs

  infTime <- get_attr(dat, "infTime") # Use infTime to find viral load in aerosol


  for (id in idsAgents) {
    if (timestep = 1) { # pseudocode
      viraLoad <- sample(possible_viral_loads, size = 1, replace = FALSE) # Pull viral load from distribution, TBD?
      dat$attributes$viralLoad[dat$agents$id == person] <- viralLoad # Assign the viral load to the person
    } else {
      if (dat$status[dat$id == id] %in% c("e", "i")) {
        dat$attributes$viralLoad[dat$agents$id == person] <- viralLoad + increment # What to increment by?
      }
    }
  }

  infectID <- cbind(active, status, infTime, viralLoad) #Just puts these attribute columns together


  return(dat)
}