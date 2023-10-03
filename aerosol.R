library(pracma)
library(EnvStats)
library(raster)
library(truncnorm)

#Set random seed for reproducibility
set.seed(19031938)
#########################
#Define Fixed Parameters

nindex = 10000

# Mask efficacy
maskfac_out = 0.0
maskfac_in = 0.0

# Minimum ach
achlimit = 0.0

# Minimum volflow per person in m3/s/per-person
volflowpplimit = 0.0

# Number of days
ndays = 12
daynumber = 1:ndays

# Number of hourly events per day
nhours = 1

# Height of room (m)
roomheight = 3.0

# Decay rate constant (1/hr); use placeholder of 1/hour
kdecay = 1.0

###################################################################################

# Integrate from 0.1 micron to 100 microns
nbins = 30
nleft = -1
nright = 2
dlogd = (nright-nleft)/nbins

dfull = numeric(nbins)

dedge = logspace(nleft, nright, nbins+1)

for (i in 1:nbins) {
  dfull[i] = (dedge[i+1]*dedge[i])^(0.5)
}

# Cappa 2 mode model for breathing - note cn here is particles/s
# Note that cn here is from Cappa 3/10/2022 4:02 am email
nmodesb  = 2

cnb = c(1.33/(1.33+0.88), 0.88/(1.33+0.88))
cmdb = c(0.45, 0.87)
gsdb = c(1.58,1.47)

# Cappa 3 mode model for talking - note cn here is particles/s
# Note that cn here is from Cappa 2/27/2022 2:52 pm email
nmodest  = 3
cnt = c(4.3/(4.3+0.1+0.011), 0.1/(4.3+0.1+0.011), 0.011/(4.3+0.1+0.011))
cmdt = c(0.45, 3.3, 10.5)
gsdt = c(1.96,1.6, 1.6)

# Set max diameter to 25 um dry ie 100 um wet
# ibinmax = 24
# Set max diameter to 5 um dry ie 20 um wet
ibinmax = 17
d = dfull[0:ibinmax]
nemb <-  array(c(numeric(ibinmax),numeric(nmodesb)),dim = c(ibinmax, nmodesb))
vemb <-  array(c(numeric(ibinmax),numeric(nmodesb)),dim = c(ibinmax, nmodesb))
nemt <-  array(c(numeric(ibinmax),numeric(nmodest)),dim = c(ibinmax, nmodest))
vemt <- array(c(numeric(ibinmax),numeric(nmodest)),dim = c(ibinmax, nmodest))



# Calculate number (1/s) and volume (um3/s) in each size bin for breathing for 1 p/s
for (im in 1:nmodesb) {
  for (i in 1:ibinmax) {
    nemb[i,im] = dlogd * log(10.0) * cnb[im]/((2.0*pi)^(0.5) * log(gsdb[im])) * exp(-((log(d[i]) - log(cmdb[im]))^2)/(2.0*(log(gsdb[im]))^2))
    vemb[i,im] = nemb[i,im] * (pi * (d[i])^3)/6.0
  }
}

# Calculate number (1/s) and volume (um3/s) in each size bin for talking for 1 p/s

for (im in 1:nmodest){
  for (i in 1:ibinmax) {
  nemt[i,im] = dlogd * log(10.0) * cnt[im]/((2.0*pi)^(0.5) * log(gsdt[im])) * exp(-((log(d[i]) - log(cmdt[im]))^2)/(2.0*(log(gsdt[im]))^2))
vemt[i,im] = nemt[i,im] * (pi * (d[i])^3)/6.0
  }
}

# Total wet volume emission rate for each bin for breathing for p/s = 1.0 in mL/s
# Scale by 64 for dry to wet volume
# Multiply by 1e-12 to convert um3/s to mL/s
vembreathe = 64.0 * 1.0e-12 * (vemb[,1]+vemb[,2])

vemtalk = 64.0 * 1.0e-12 * (vemt[,1]+vemt[,2]+vemt[,3])

# Specify distribution of p/s for breathing
# From Cappa's 3/11 and 3/12 emails (also see my email to him on 3/17 at 2.40 pm) for breathing
# Scale by factor of 8 (ie add ln(8) to mu) to account for APS sampling
mub = log(8) - 0.48*log(10)
sigmab = 0.85*log(10)/(2^(0.5))
nemm_breathe <- scan("data/data3.csv", sep=",")


# Specify distribution of p/s for talking
# From Asadi, fig 4B, draw samples of number of particles/s emitted for talking
# Mean = 3.4 particles/s, Std = 2.7 particles/s
# Scale by factor of 13 (ie multiply mean and std dev by 13) to account for APS sampling
# Truncate log-normal at 260 particles/s (limit of 20 from Fig 5A in Asadi et al. scaled by factor of 13)

mean = 3.4*13
std = 2.7*13
mu = log(mean*mean/((mean*mean + std*std)^(0.5)))
sigma = (log(1.0 + std*std/(mean*mean)))^0.5
lmax = log(20.0*13)
bclip = (lmax - mu)/sigma
lclip <- scan("data/data2.csv", sep=",")
nemm_talk = exp(lclip)

#lclip = truncnorm.rvs(-np.inf, bclip, loc=mu, scale=sigma, size=nindex,random_state=rng)

# Set p/s from breathing to no more than 10% from talking
for (i in 1:nindex) {
  if (nemm_breathe[i] > 0.1*nemm_talk[i]) {
    nemm_breathe[i] = 0.1*nemm_talk[i]
  }
}

# Set p/s from breathing to a max of 24 p/s (limit of 3 from Fig 5A in Asadi et al. scaled by factor of 8)
nemm_breathe = clamp(nemm_breathe, lower = -Inf, 24, useValue = TRUE)
sprintf('nemm_talk, %f, %f, %f, %f, %f', mean(nemm_talk), median(nemm_talk), std(nemm_talk), min(nemm_talk), max(nemm_talk))
sprintf('nemm_breathe, %f, %f, %f, %f, %f', mean(nemm_breathe),median(nemm_breathe),std(nemm_breathe),min(nemm_breathe),max(nemm_breathe))
#######################################

vbreatheindex = numeric(nindex)
vtalkindex = numeric(nindex)

for (n in 1:nindex) {
  vbreatheindex[n] = nemm_breathe[n] * sum(vembreathe)
  vtalkindex[n] = nemm_talk[n] * sum(vemtalk)
}

vbreatheindex_median = median(vbreatheindex)
vtalkindex_median = median(vtalkindex)
vbreatheindex_mean = mean(vbreatheindex)
vtalkindex_mean = mean(vtalkindex)

vtalkindex_75 = quantile(vtalkindex, probs=.75)
vtalkindex_90 = quantile(vtalkindex, probs=.90)

vratio = vtalkindex/vbreatheindex

############
# BREATHING RATE

# Breathing rate in m3/s
# Uniform distribution from sitting to light activity from Table 3 Henriques et al. (2022)

# Breathing rate in m3/s
# Uniform distribution from sitting to light activity from Table 3 Henriques et al. (2022)

vdotbreathe = array(c(numeric(nindex),numeric(ndays)),dim = c(nindex, ndays))

lowbreathe = 0.51/3600.0
highbreathe = 1.24/3600.0
# lowbreathe = 1.0e-4
# highbreathe = 2.0e-4

br_activity = c(0.51/3600.0, 0.57/3600.0, 1.24/3600)
activityindex = array(c(numeric(nindex), numeric(ndays)), dim = c(nindex, ndays))

for (iday in 1:ndays) {
  activityindex[,iday] = sample(3, size = nindex, replace = TRUE)

  for (iindex in 1:nindex) {
    vdotbreathe[iindex,iday]=br_activity[activityindex[iindex,iday]]
  }
}


############
# VIRAL LOAD AND INFECTIOUSNESS

# Viral load NP (copies)

vnose <- matrix(0, nrow = nindex, ncol = ndays)
kevn <- read.csv("data/Ke_Log10VLs_10000.csv", sep=",", skip = 1)

for (iday in 1:ndays) {
  irow <- 10 * iday + 5
  vnose[, iday] <- unlist((10^kevn[irow, 2:10001]) * 3)
}


vsalv = matrix(0, nrow = nindex, ncol = ndays)
kevs = read.csv("data/Ke_Log10VLs_saliva_10000.csv", sep=",", skip = 1)
for (iday in 1:ndays) {
  irow <- 10 * iday + 5
  vnose[, iday] <- unlist((10^kevs[irow, 2:10001]))
}


vinfectiousness = matrix(0, nrow = nindex, ncol = ndays)
kevi = read.csv("data/Ke_Infectiousness_10000.csv", sep=",", skip = 1)
for (iday in 1:ndays) {
  irow <- 10 * iday + 5
  vnose[, iday] <- unlist((3*kevi[irow, 2:10001]))
}


##############
# ACH
ach = matrix(0, nrow = nindex, ncol = ndays)
achm = matrix(0, nrow = nindex, ncol = ndays)

# From Persily et al. (2005)
meanach = 2.00
stdach = 2.45

mu_ach = log(meanach*meanach/(meanach*meanach + stdach*stdach)^0.5)
sigma_ach = (log(1.0 + stdach*stdach/(meanach*meanach)))^0.5

for (iday in 1:ndays){ ###########  May need to replace this with saved data for replicability, or rely only on R seed
  ach[,iday] = rlnorm(mu_ach,sigma_ach, n=nindex)
  achm[, iday] = clamp(ach[, iday], achlimit, Inf)
}






