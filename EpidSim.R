### DKP note
## see https://statnet.org/nme/ and https://statnet.org/nme/d3.html and https://statnet.org/nme/d4-s9-SEIR.html

## Code based on:
##
## Network Modeling for Epidemics 2015
## Day 5 | Tutorial 3: An SEIR Epidemic
##



rm(list = ls())

#install.packages("EpiModel")
library("EpiModel")


set.seed(27708)


### call of specify the network here if you'd like to be more specific
### (otherwise, set with basic preferential age-mixing below)


# modified to add an E for SEIR -- eventually should be SEIRS depending on sim length
### Transmission module / adding an E compartment





infect <- function(dat, at) {

  ## Uncomment this to run environment interactively
  #browser()

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")

  ## Parameters ##
  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")

  ## Find infected nodes ##
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  ## Initialize default incidence at 0 ##
  nInf <- 0

  ## If any infected nodes, proceed with transmission ##
  if (nElig > 0 && nElig < nActive) {

    ## Look up discordant edgelist ##
    del <- discord_edgelist(dat, at)

    ## If any discordant pairs, proceed ##
    if (!(is.null(del))) {

      # Set parameters on discordant edgelist data frame
      del$transProb <- inf.prob
      del$actRate <- act.rate
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate

      # Stochastic transmission process
      transmit <- rbinom(nrow(del), 1, del$finalProb)

      # Keep rows where transmission occurred
      del <- del[which(transmit == 1), ]

      # Look up new ids if any transmissions occurred
      idsNewInf <- unique(del$sus)
      nInf <- length(idsNewInf)

      # Set new attributes and transmission matrix
      if (nInf > 0) {
        status[idsNewInf] <- "e"
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "status", status)
        dat <- set_attr(dat, "infTime", infTime)

        dat <- set_transmat(dat, del, at)

      }

      # Set move from recovery to susceptibility after 8 months
      if (nInf > 0) {
        status[idsNewInf] <- "e"
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "status", status)
        dat <- set_attr(dat, "infTime", infTime)

        dat <- set_transmat(dat, del, at)

      }

    }
  }

  ## Save summary statistic for S->E flow
  dat <- set_epi(dat, "se.flow", at, nInf)

  ## Save summary statistic for R->S flow
  dat <- set_epi(dat, "rs.flow", at, nInf)

  return(dat)
}


progress <- function(dat, at) {   # also with an E compartment

  ## Uncomment this to function environment interactively
  #browser()

  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

  ## Parameters ##
  ei.rate <- get_param(dat, "ei.rate")
  ir.rate <- get_param(dat, "ir.rate")
  rs.rate <- get_param(dat, "rs.rate")  # needed for SEIRS

  ## E to I progression process ##
  nInf <- 0
  idsEligInf <- which(active == 1 & status == "e")
  nEligInf <- length(idsEligInf)

  if (nEligInf > 0) {
    vecInf <- which(rbinom(nEligInf, 1, ei.rate) == 1)
    if (length(vecInf) > 0) {
      idsInf <- idsEligInf[vecInf]
      nInf <- length(idsInf)
      status[idsInf] <- "i"
    }
  }

  ## I to R progression process ##
  nRec <- 0
  idsEligRec <- which(active == 1 & status == "i")
  nEligRec <- length(idsEligRec)

  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, ir.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nRec <- length(idsRec)
      status[idsRec] <- "r"
    }
  }


  ## R to S progression process ##
  nSus <- 0
  idsEligSus <- which(active == 1 & status == "r")
  nEligSus <- length(idsEligSus)

  if (nEligSus > 0) {
    vecSus <- which(rbinom(nEligSus, 1, rs.rate) == 1)
    if (length(vecSus) > 0) {
      idsSus <- idsEligRec[vecSus]
      nSus <- length(idsSus)
      status[idsSus] <- "s"
    }
  }

  ## Write out updated status attribute ##
  dat <- set_attr(dat, "status", status)

  ## Save summary statistics ##
  dat <- set_epi(dat, "ei.flow", at, nInf)
  dat <- set_epi(dat, "ir.flow", at, nRec)
  dat <- set_epi(dat, "rs.flow", at, nSus)
  dat <- set_epi(dat, "e.num", at,
                 sum(active == 1 & status == "e"))
  dat <- set_epi(dat, "r.num", at,
                 sum(active == 1 & status == "r"))
  dat <- set_epi(dat, "s.num", at,
                 sum(active == 1 & status == "s"))

  return(dat)
}



xpinf_1 <- read.csv("data/xpinf_1.csv")


# Create a `list.of.updaters`
#Will replace this with a function that builds the list (with num tsteps as var)
list.of.updaters <- list(
  # this is one updater
  list(
    at = 1,
    param = list(
      inf.prob = sample(xpinf_1[, 1], size = 1)
    )
  ),
  # this is another updater
  list(
    at = 2,
    param = list(
      inf.prob = sample(xpinf_1[, 2], size = 1)
    )
  ),
  list(
    at = 3,
    param = list(
      inf.prob = sample(xpinf_1[, 3], size = 1)
    )
  ),
  list(
    at = 4,
    param = list(
      inf.prob = sample(xpinf_1[, 4], size = 1)
    )
  ),
  list(
    at = 5,
    param = list(
      inf.prob = sample(xpinf_1[, 5], size = 1)
    )
  ),
  list(
    at = 6,
    param = list(
      inf.prob = sample(xpinf_1[, 6], size = 1)
    )
  ),
  list(
    at = 7,
    param = list(
      inf.prob = sample(xpinf_1[, 7], size = 1)
    )
  ),
  list(
    at = 8,
    param = list(
      inf.prob = sample(xpinf_1[, 8], size = 1)
    )
  ),
  list(
    at = 9,
    param = list(
      inf.prob = sample(xpinf_1[, 9], size = 1)
    )
  ),
  list(
    at = 10,
    param = list(
      inf.prob = sample(xpinf_1[, 10], size = 1)
    )
  ),
  list(
    at = 11,
    param = list(
      inf.prob = sample(xpinf_1[, 11], size = 1)
    )
  ),
  list(
    at = 12,
    param = list(
      inf.prob = sample(xpinf_1[, 12], size = 1)
    )
  )
)


param <- param.net(inf.prob = 0.5, act.rate = 2, ei.rate = 0.01, ir.rate = 0.01,
                   rs.rate = 0.01, .param.updater.list = list.of.updaters)

init <- init.net(i.num = 10, status.rand = FALSE)

# Ignore the type = "SI" setting; this is a bug that will be fixed
control <- control.net(type = NULL, nsteps = 12, nsims = 1,
                       infection.FUN = infect, progress.FUN = progress,
                       recovery.FUN = NULL, skip.check = TRUE,
                       resimulate.network = FALSE, verbose.int = 0,
                       verbose = TRUE,
                       save.transmat = TRUE, save.other = c("attr"),
                       summary_nets.FUN = summary_nets)

##WE WANT TO WRITE NEW infection.fun
#The nodal attributes will eventually be utilized to weight the distribution into rooms
#But as of now, we are just going to focus on getting into rooms and applying
#aerosol



## Network model
nw <- network.initialize(100, directed = FALSE)  ## this is a 500-person network with no preferential attachment
# set age-mixing:  https://rpubs.com/smjenness/111058
#nw <- network.initialize(500, directed = FALSE, )  ## this uses the network specified above
age <- sample(1:99, 1000, TRUE)

agecat <- rep(2, 1000)
agecat[age < 18] <- 1
agecat[age >= 65] <- 3
table(agecat)
nw <- set.vertex.attribute(nw, "age", age)
nw <- set.vertex.attribute(nw, "agecat", agecat)
nw

mean.deg <- 0.75
edges <- mean.deg * (1000/2)


# adds age-mixing to network
nodematch.targets <- c(120, 80, 100)  # specifies most age-clustering in younger ages, less in middle, more in upper
fit.agecat <- ergm(nw ~ edges + nodematch("agecat", diff = TRUE),
                   target.stats = c(edges, nodematch.targets))
sim.agecat <- simulate(fit.agecat, nsim = 1e4, statsonly = TRUE)
colMeans(sim.agecat)
# plot below shows clustering by age
el <- as.edgelist(simulate(fit.agecat))
age.el <- matrix(age[el], ncol = 2)
plot(jitter(age.el[, 1]), jitter(age.el[, 2]))
abline(a = 0, b = 1)

# estimate network parameters
est <- netest(nw, formation = ~edges, target.stats = 150,
              coef.diss = dissolution_coefs(~offset(edges), 10))  # this is the generic network
est <- netest(nw, formation = ~edges + nodematch("agecat", diff = TRUE),
              target.stats = c(edges, nodematch.targets),
              coef.diss = dissolution_coefs(~offset(edges), 10))  # this is the network with age-mixing
est <- netest(nw, formation1, target.stats1, coef.diss,
              coef.form = -Inf, set.control.ergm = control.ergm(MCMLE.maxit = 500))


## Simulate
sim <- netsim(est, param, init, control)

## Plot
par(mar = c(3,3,1,1), mgp = c(2,1,0))
#plot(sim, y = c("s.num", "i.num", "e.num", "r.num"), popfrac = FALSE,
#     mean.col = 1:4, qnts = 1, qnts.col = 1:4, leg = TRUE)

# ^^ This plot gave errors ^^

plot(sim)
plot(sim, y = c("se.flow", "ei.flow", "ir.flow","rs.flow"), legend = TRUE)

## next:  extend to SEIRS:  https://www.nature.com/articles/s41592-020-0856-2

