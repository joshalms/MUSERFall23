assignRooms <- function(idsAgents) {

  nRooms = 10

  Rooms <- vector("list", length = nRooms)

  for (id in idsAgents) {
    idx <- sample(1:nRooms, 1)
    Rooms[[idx]] <- c(Rooms[[idx]], id)
  }

  return(Rooms)

}