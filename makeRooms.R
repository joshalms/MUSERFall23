assignRooms <- function(idsAgents) {

  nRooms = 10 # Number of rooms

  Rooms <- vector("list", length = nRooms) # Creates list of specified number of rooms

  for (id in idsAgents) { # Iterate through each person using their ID
    idx <- sample(1:nRooms, 1) # Chooses a random room to add person to
    Rooms[[idx]] <- c(Rooms[[idx]], id) # Adds that person's ID to the room
  }

  return(Rooms) # Returns updated list of rooms. Each room is represented as a list containing the IDs of the people assigned to that room

}
