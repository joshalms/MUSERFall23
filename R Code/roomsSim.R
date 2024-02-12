

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

  source("getRoomProb.R") #Load in the function in this file
  #Right now, it is just a static inf.prob. Want to replace with output of
  #aerosol model
  source("makeRooms.R") #Load function. Splits all ids into 10 rooms (10 is arbitrary constant that we can change)s



  ## Attributes ##
  active <- get_attr(dat, "active") #list of all people who are doing something at this time step
  #Idea is that some ppl won't interact with anyone
  status <- get_attr(dat, "status") #Lit of all statuses with indices corresponding to ids

  infTime <- get_attr(dat, "infTime") # Use infTime to find viral load in aerosol

  nodes <- cbind(active, status, infTime) #Just puts these attribute columns together

  ## Find infected nodes ##

  idsAgents <- which(active == 1) #Get ids for all agents who are doing something
  nActive <- sum(active == 1)

  ## Initialize default incidence at 0 ##
  nInf <- 0

# FOR:
    #Assign each person to k rooms (we define k)

    #Call the aerosol model on each of these rooms
      #Add a vertex attribute for viral load

    #Modify the aerosol script to be a function

      #Run as usual but take a list of people with viral loads and output
      #One probability - sample
      #Sample from a binomial distribution approximated w/ mean of that probability
      #Randomly infect that number of people in the room


  ## If any infected nodes, proceed with transmission ##

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

    Rooms <- assignRooms(idsAgents) #Split people into 10 rooms

    for (room in Rooms) { #Iterate through the rooms

      # ** Load kevn, kevs, kevi

      infProb <- getRoomProb() #Need to restructure aerosol to give inf prob here
      sus <- room[which(status == "s")] #Create a vector that holds all person ids for the susceptible in the room
      transmissions <- rbinom(length(sus), 1, infProb) #Sample from binomial distribution for number of infections
      idsNewInf <- sus[which(transmissions == 1)] #New vector that is a subset of sus -> contains only ids for infected
      nInf <- length(idsNewInf) #Number of infected for the room


      if (nInf > 0) { #If anyone was infected

        toInfect <- c()

        for (person in idsNewInf) {
            toInfect <- c(toInfect, person) #Makes a new vector of the ids we want
        }

               status[toInfect] <- "e" #Update to exposed
               infTime[toInfect] <- at #Set infection time to current time
               dat <- set_attr(dat, "status", status)
               dat <- set_attr(dat, "infTime", infTime) #Update model state
             }

    }




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
   dat <- set_epi(dat, "se.flow", at, nInf)
  #
  # ## Save summary statistic for R->S flow
   dat <- set_epi(dat, "rs.flow", at, nInf)




  return(dat)












  return(dat)

}
