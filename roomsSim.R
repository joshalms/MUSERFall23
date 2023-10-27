

#This is the primary function to simulate infections and transmission of
#an aerosol-based disease

#Approach:
# At each time step we randomly assign every person to exist within k spaces
# out of n total for k < n.
# Within each of these simulated bins, apply the aerosol model based on the
# number of people in the room and the viral load of Is present
# Based on output of aerosol for x infected, randomly assign x of the agents
# within the room to transition states into "E"



infect <- function(dat, at) {


  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")

  ## Find infected nodes ##
  idsInf <- which(active == 1 & status == "i")
  idsAgents <- which(active == 1) #Get ids for all agents
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  ## Initialize default incidence at 0 ##
  nInf <- 0




  ## If any infected nodes, proceed with transmission ##
  if (nElig > 0 && nElig < nActive) {

    #Assign agents to rooms
    #Call aerosol model on each room and get random variable number of
    #transmissions per room
    #increment nInf
    #Randomly assign x s agents in a room to transition -> repeat for each room

    #For agent in idsAgents length:
      #Randomly append to k of n rooms (list of lists)

    #Iterate through rooms:
      #Call aerosol model
      #Randomly assign transitions

    #Update status & params

    ## Look up discordant edgelist ##
  #   del <- discord_edgelist(dat, at)
  #
  #   ## If any discordant pairs, proceed ##
  #   if (!(is.null(del))) {
  #
  #     # Set parameters on discordant edgelist data frame
  #     del$transProb <- inf.prob
  #     del$actRate <- act.rate
  #     del$finalProb <- 1 - (1 - del$transProb)^del$actRate
  #
  #     # Stochastic transmission process
  #     transmit <- rbinom(nrow(del), 1, del$finalProb)
  #
  #     # Keep rows where transmission occurred
  #     del <- del[which(transmit == 1), ]
  #
  #     # Look up new ids if any transmissions occurred
  #     idsNewInf <- unique(del$sus)
  #     nInf <- length(idsNewInf)
  #
  #     # Set new attributes and transmission matrix
  #     if (nInf > 0) {
  #       status[idsNewInf] <- "e"
  #       infTime[idsNewInf] <- at
  #       dat <- set_attr(dat, "status", status)
  #       dat <- set_attr(dat, "infTime", infTime)
  #
  #       dat <- set_transmat(dat, del, at)
  #
  #     }
  #
  #     # Set move from recovery to susceptibility after 8 months
  #     if (nInf > 0) {
  #       status[idsNewInf] <- "e"
  #       infTime[idsNewInf] <- at
  #       dat <- set_attr(dat, "status", status)
  #       dat <- set_attr(dat, "infTime", infTime)
  #
  #       dat <- set_transmat(dat, del, at)
  #
  #     }
  #
  #   }
  # }
  #
  # ## Save summary statistic for S->E flow
  # dat <- set_epi(dat, "se.flow", at, nInf)
  #
  # ## Save summary statistic for R->S flow
  # dat <- set_epi(dat, "rs.flow", at, nInf)

  return(dat)
}










  return(dat)

}
