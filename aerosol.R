library(pracma)
library(EnvStats)
library(raster)
library(truncnorm)
library(msm)

#Set random seed for reproducibility
set.seed(19031938)
float_epsilon <- 2.220446049250313e-16
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
nemm_breathe <- rlnorm(mub,sigmab, n=nindex) #rng.lognormal(mub, sigmab, size=nindex)
#rlnorm(mu_ach,sigma_ach, n=nindex)


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
bclip <-  (lmax - mu)/sigma
lclip <- scan("data/data2.csv", sep=",")
nemm_talk = exp(lclip)

#lclip = truncnorm.rvs(-np.inf, bclip, loc=mu, scale=sigma, size=nindex,random_state=rng)
#rtruncnorm(n = nindex, a = aclip, b = bclip, mean = mu_tend, sd = sigma_tend)


# Set p/s from breathing to no more than 10% from talking
for (i in 1:nindex) {
  if (nemm_breathe[i] > 0.1*nemm_talk[i]) {
    nemm_breathe[i] = 0.1*nemm_talk[i]
  }
}

# Set p/s from breathing to a max of 24 p/s (limit of 3 from Fig 5A in Asadi et al. scaled by factor of 8)
nemm_breathe = clamp(nemm_breathe, lower = -Inf, upper = 24, useValues = TRUE)
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


##############
# DEPOSITION

# Calculate diameter-dependent deposition velocity based on Lai and Nazaroff (2000)
#

# Air viscosity (kg/m-s); Henriques (2022)
mu_air = 1.8e-5

lamda_air = 0.0651

# Particle density (kg/m3); Henriques (2022)
dens_particle = 1000.0

# Air density (kg/m3); Henriques (2022)
dens_air = 1.2

# Gravitational constant (m/s2)
grav = 9.81
#
# Knudsen number
Kn = 2.0 * lamda_air/d

Cc = 1.0 + Kn * (1.257 + 0.4 * exp(-1.1/Kn))

setvel_particle = ((d * 1.0e-6)**2) * (dens_particle - dens_air) * grav * Cc /(18.0 * mu_air)

ustar = 0.01

# Air kinematic viscosity (m2/s)
kinem_air = mu_air/dens_air
#
# Boltzman constant (m2 kg/s2-K)
boltz = 1.38e-23

# temperature (K)
temperature = 298.0

# Particle diffusivity in air (m2/s); Equation 8.73, Seinfeld and Pandis (1998)
diff_particle = boltz * temperature * Cc /(3 * pi * mu_air * d*1.0e-6)

# Schmidt number
Sc  = kinem_air/diff_particle

# Calculate deposition velocity (m/s) for upward facing horizontal surface; Table 2, Lai and Nazaroff (2000)

# rplus; Lai and Nazaroff (2000)
rplus = d*1.0e-6/2.0 * ustar/kinem_air

scpow = Sc^(-1.0/3.0)

aconst = 0.5 * log(((10.92*scpow + 4.3)^3)/(1.0/Sc + 0.0609)) + (3)^0.5 * atan((8.6 - 10.92*scpow)/(3^0.5 * 10.92 * scpow))

bconst = 0.5 * log(((10.92*scpow + rplus)^3)/(1.0/Sc + 7.669e-4*(rplus**3))) + 3^0.5 * atan((2.0*rplus - 10.92*scpow)/(3^0.5 * 10.92 * scpow))

integ = 3.64 * (Sc**0.66667) * (aconst - bconst) + 39.0

# Deposition velocity (m/s)
dep = setvel_particle/(1.0 - exp(-setvel_particle*integ/ustar))


# Deposition rate coefficient (1/hr)
kdep = 3600.0*dep/roomheight

##############
# DIAMETER-DEPENDENT DEPOSITION FRACTION IN RESPIRATORY SYSTEM - from Henriques et al. (2022)

ifrac = 1.0 - 0.5*(1.0 - 1.0/(1.0 + 0.00076*(d^2.8)))
fdep = ifrac * (0.0587 + 0.911/(1.0 + exp(4.77 + 1.485*log(d))) + 0.943/(1.0 + exp(0.508 - 2.58*log(d))))

##############

# AREA PER-PERSON (m**2/person)

areapp = matrix(0, nrow = nindex, ncol = ndays)
densx = matrix(0, nrow = nindex, ncol = ndays)

# ASHRAE default densities

# Generate random values for areapp
for (iday in 1:ndays) {
  areapp[, iday] <- runif(nindex, min = 1.0, max = 10.0)
}

# Calculate densx values
for (iday in 1:ndays) { #May again need to replace this with loaded data
  densx[, iday] <- 100.0 / areapp[, iday]
}


######
# NUMBER OF SUSCEPTIBLE CONTACTS

nsus = matrix(0, nrow = nindex, ncol = ndays)
for (iday in 1:ndays) {
  nsus[, iday] <- round(rgamma(n = nindex, shape = 20.0/30.0, scale = 30))
}
nsus = apply(nsus, MARGIN = c(1,2), FUN = as.integer)
# changed (AB) - previously converted entire matrix into some int list?


###############
# EVENT DURATION (s)

tend = matrix(0, nrow = nindex, ncol = ndays)


mu_tend = log(1.0)


sigma_tend=0.5
lntend_min = log(0.1)
lntend_max = log(3.0+float_epsilon)

aclip = (lntend_min - mu_tend)/sigma_tend
bclip = (lntend_max - mu_tend)/sigma_tend

#Random element that may need to be replaced with the corresponding file
for (iday in 1:ndays) {
  lntend_clip <- rtruncnorm(n = nindex, a = aclip, b = bclip, mean = mu_tend, sd = sigma_tend)
  tend[, iday] <- 3600.0 * exp(lntend_clip)
}

###############
# TALK TIME FRACTION

talkfraction = matrix(0, nrow = nindex, ncol = ndays)

for (iday in 1:ndays) {
  talkfraction[, iday] <- runif(nindex, min = 0.05, max = 0.25 + float_epsilon)
}

##############
# CALCULATE AREA OF ROOM (m**2)
roomarea = (1 + nsus)*areapp

# CALCULATE VOLUME OF ROOM (m**3)

roomvol = roomheight*roomarea

###############
# CALCULATE VOLUME FLOW RATE (m3/s)

volflow = roomvol * ach / 3600.0
volflowmin = 1e-3 * (5.0 * (1 + nsus) + 1.0 * roomarea)
volflownew = pmax(volflow, volflowmin) # clamp gave error for me (AB)
achm = volflownew * 3600.0 / roomvol
volflowpp = volflow / (1 + nsus)
volflowppm = volflownew / (1 + nsus)

###############
# MAIN CALCULATIONS
# const_1, xconst_1, const_2, xconst_2, const_3, xconst_3 chosen by trial and error to get R0 = 3
# const_1, xconst_1 still to be adjusted

# Linear - SALIVA  - ONLY FAR FIELD INTERACTIONS
hpower_1 = 1.0
const_1 = 0.4498
# Linear - SALIVA  - NEAR+FAR FIELD INTERACTIONS
xconst_1 = 0.216
xdil_1 = 24

# Linear - NOSE - ONLY FAR FIELD INTERACTIONS
hpower_2 = hpower_1
const_2 = 2.049
# Linear - NOSE - NEAR+FAR FIELD INTERACTIONS
xconst_2 = 0.9677
xdil_2 = 24

# Power law - SALIVA - ONLY FAR FIELD INTERACTIONS
hpower_3 = runif(nindex, min = 0.5, max = 0.75)
const_3 = 271.90
# Power law - SALIVA - NEAR+FAR FIELD INTERACTIONS
xconst_3 = 128.28
xdil_3 = 24

# CAUTION - const_4 and xconst_4 have not yet been adjusted to get R0 = 3
# Power law 0.6 - NOSE - ONLY FAR FIELD INTERACTIONS
hpower_4 = hpower_3
const_4 = 2.673e3
# Power law 0.6 - NOSE - NEAR+FAR FIELD INTERACTIONS
xconst_4 = 1.3355e3
xdil_4 = 24

# Infection virions per mL of aerosol volume
# Note vsalv is per ml of saliva, vnose is per mL of VMT, vinfectioiusness is per swab
vinf_1 = const_1 * (vsalv ** hpower_1)
vinf_2 = const_2 * (vnose ** hpower_2)
vinf_3 = matrix(0, nrow = nindex, ncol = ndays)
vinf_4 = matrix(0, nrow = nindex, ncol = ndays)
for (iday in 1:ndays) {
  vinf_3[, iday] <- const_3 * (vsalv[, iday] ^ hpower_3)
  vinf_4[, iday] <- const_4 * (vnose[, iday] ^ hpower_4)
}

xvinf_1 = xconst_1 * (vsalv ^ hpower_1)
xvinf_2 = xconst_2 * (vnose ^ hpower_2)
xvinf_3 = matrix(0, nrow = nindex, ncol = ndays)
xvinf_4 = matrix(0, nrow = nindex, ncol = ndays)
for (iday in 1:ndays) {
  xvinf_3[, iday] <- xconst_3 * (vsalv[, iday] ^ hpower_3)
  xvinf_4[, iday] <- xconst_4 * (vnose[, iday] ^ hpower_4)
}

# ID50 in virions
# SET TO 1 to INCLUDE IN SCALING CONSTANT
id50 = 1.0

# ARRAYS FOR WELL-MIXED CASE
rnaconc_1 = matrix(0, nrow = ibinmax, ncol = ndays)
rnaconc_2 = matrix(0, nrow = ibinmax, ncol = ndays)
rnaconc_3 = matrix(0, nrow = ibinmax, ncol = ndays)
rnaconc_4 = matrix(0, nrow = ibinmax, ncol = ndays)
dose_1 = matrix(0, nrow = nindex, ncol = ndays)
dose_2 = matrix(0, nrow = nindex, ncol = ndays)
dose_3 = matrix(0, nrow = nindex, ncol = ndays)
dose_4 = matrix(0, nrow = nindex, ncol = ndays)
pinf_1 = matrix(0, nrow = nindex, ncol = ndays)
pinf_2 = matrix(0, nrow = nindex, ncol = ndays)
pinf_3 = matrix(0, nrow = nindex, ncol = ndays)
pinf_4 = matrix(0, nrow = nindex, ncol = ndays)
ninf_1 = matrix(0, nrow = nindex, ncol = ndays)
ninf_2 = matrix(0, nrow = nindex, ncol = ndays)
ninf_3 = matrix(0, nrow = nindex, ncol = ndays)
ninf_4 = matrix(0, nrow = nindex, ncol = ndays)
ninf_1x = matrix(0, nrow = nindex, ncol = ndays)
ninf_2x = matrix(0, nrow = nindex, ncol = ndays)
ninf_3x = matrix(0, nrow = nindex, ncol = ndays)
ninf_4x = matrix(0, nrow = nindex, ncol = ndays)
Hr = matrix(0, nrow = nindex, ncol = ndays)
HrV = matrix(0, nrow = nindex, ncol = ndays)
Hrtrue = matrix(0, nrow = nindex, ncol = ndays)
vbreathetrue = matrix(0, nrow = nindex, ncol = ndays)
vtalktrue = matrix(0, nrow = nindex, ncol = ndays)

rnaconcm_1 = matrix(0, nrow = ibinmax, ncol = ndays)
rnaconcm_2 = matrix(0, nrow = ibinmax, ncol = ndays)
rnaconcm_3 = matrix(0, nrow = ibinmax, ncol = ndays)
rnaconcm_4 = matrix(0, nrow = ibinmax, ncol = ndays)
dosem_1 = matrix(0, nrow = nindex, ncol = ndays)
dosem_2 = matrix(0, nrow = nindex, ncol = ndays)
dosem_3 = matrix(0, nrow = nindex, ncol = ndays)
dosem_4 = matrix(0, nrow = nindex, ncol = ndays)
pinfm_1 = matrix(0, nrow = nindex, ncol = ndays)
pinfm_2 = matrix(0, nrow = nindex, ncol = ndays)
pinfm_3 = matrix(0, nrow = nindex, ncol = ndays)
pinfm_4 = matrix(0, nrow = nindex, ncol = ndays)
ninfm_1 = matrix(0, nrow = nindex, ncol = ndays)
ninfm_2 = matrix(0, nrow = nindex, ncol = ndays)
ninfm_3 = matrix(0, nrow = nindex, ncol = ndays)
ninfm_4 = matrix(0, nrow = nindex, ncol = ndays)
ninfm_1x = matrix(0, nrow = nindex, ncol = ndays)
ninfm_2x = matrix(0, nrow = nindex, ncol = ndays)
ninfm_3x = matrix(0, nrow = nindex, ncol = ndays)
ninfm_4x = matrix(0, nrow = nindex, ncol = ndays)

# ARRAYS FOR NEAR+FAR CASE

nsus_near = matrix(0, nrow = nindex, ncol = ndays)
fnear = matrix(0, nrow = nindex, ncol = ndays)
tnear = matrix(0, nrow = nindex, ncol = ndays)

xrnaconc_far_1 = array(0, dim = c(ibinmax, nindex, ndays))
xrnaconc_far_2 = array(0, dim = c(ibinmax, nindex, ndays))
xrnaconc_far_3 = array(0, dim = c(ibinmax, nindex, ndays))
xrnaconc_far_4 = array(0, dim = c(ibinmax, nindex, ndays))
xrnaconc_near_1 = array(0, dim = c(ibinmax, nindex, ndays))
xrnaconc_near_2 = array(0, dim = c(ibinmax, nindex, ndays))
xrnaconc_near_3 = array(0, dim = c(ibinmax, nindex, ndays))
xrnaconc_near_4 = array(0, dim = c(ibinmax, nindex, ndays))
xdose_far_1 = matrix(0, nrow = nindex, ncol = ndays)
xdose_far_2 = matrix(0, nrow = nindex, ncol = ndays)
xdose_far_3 = matrix(0, nrow = nindex, ncol = ndays)
xdose_far_4 = matrix(0, nrow = nindex, ncol = ndays)
xdose_nearfar_1 = matrix(0, nrow = nindex, ncol = ndays)
xdose_nearfar_2 = matrix(0, nrow = nindex, ncol = ndays)
xdose_nearfar_3 = matrix(0, nrow = nindex, ncol = ndays)
xdose_nearfar_4 = matrix(0, nrow = nindex, ncol = ndays)
xpinf_far_1 = matrix(0, nrow = nindex, ncol = ndays)
xpinf_far_2 = matrix(0, nrow = nindex, ncol = ndays)
xpinf_far_3 = matrix(0, nrow = nindex, ncol = ndays)
xpinf_far_4 = matrix(0, nrow = nindex, ncol = ndays)
xpinf_nearfar_1 = matrix(0, nrow = nindex, ncol = ndays)
xpinf_nearfar_2 = matrix(0, nrow = nindex, ncol = ndays)
xpinf_nearfar_3 = matrix(0, nrow = nindex, ncol = ndays)
xpinf_nearfar_4 = matrix(0, nrow = nindex, ncol = ndays)
xpinf_1 = matrix(0, nrow = nindex, ncol = ndays)
xpinf_2 = matrix(0, nrow = nindex, ncol = ndays)
xpinf_3 = matrix(0, nrow = nindex, ncol = ndays)
xpinf_4 = matrix(0, nrow = nindex, ncol = ndays)
xninf_1 = matrix(0, nrow = nindex, ncol = ndays)
xninf_2 = matrix(0, nrow = nindex, ncol = ndays)
xninf_3 = matrix(0, nrow = nindex, ncol = ndays)
xninf_4 = matrix(0, nrow = nindex, ncol = ndays)

xrnaconcm_far_1 = array(0, dim = c(ibinmax, nindex, ndays))
xrnaconcm_far_2 = array(0, dim = c(ibinmax, nindex, ndays))
xrnaconcm_far_3 = array(0, dim = c(ibinmax, nindex, ndays))
xrnaconcm_far_4 = array(0, dim = c(ibinmax, nindex, ndays))
xrnaconcm_near_1 = array(0, dim = c(ibinmax, nindex, ndays))
xrnaconcm_near_2 = array(0, dim = c(ibinmax, nindex, ndays))
xrnaconcm_near_3 = array(0, dim = c(ibinmax, nindex, ndays))
xrnaconcm_near_4 = array(0, dim = c(ibinmax, nindex, ndays))
xdosem_far_1 = matrix(0, nrow = nindex, ncol = ndays)
xdosem_far_2 = matrix(0, nrow = nindex, ncol = ndays)
xdosem_far_3 = matrix(0, nrow = nindex, ncol = ndays)
xdosem_far_4 = matrix(0, nrow = nindex, ncol = ndays)
xdosem_nearfar_1 = matrix(0, nrow = nindex, ncol = ndays)
xdosem_nearfar_2 = matrix(0, nrow = nindex, ncol = ndays)
xdosem_nearfar_3 = matrix(0, nrow = nindex, ncol = ndays)
xdosem_nearfar_4 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_far_1 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_far_2 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_far_3 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_far_4 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_nearfar_1 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_nearfar_2 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_nearfar_3 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_nearfar_4 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_far_1 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_far_2 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_far_3 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_far_4 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_nearfar_1 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_nearfar_2 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_nearfar_3 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_nearfar_4 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_1 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_2 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_3 = matrix(0, nrow = nindex, ncol = ndays)
xpinfm_4 = matrix(0, nrow = nindex, ncol = ndays)
xninfm_1 = matrix(0, nrow = nindex, ncol = ndays)
xninfm_2 = matrix(0, nrow = nindex, ncol = ndays)
xninfm_3 = matrix(0, nrow = nindex, ncol = ndays)
xninfm_4 = matrix(0, nrow = nindex, ncol = ndays)

Hrm = matrix(0, nrow = nindex, ncol = ndays)

vaerosoltot = matrix(0, nrow = nindex, ncol = ndays)
vaerosoltottalk = matrix(0, nrow = nindex, ncol = ndays)
vaerosoltotbreathe = matrix(0, nrow = nindex, ncol = ndays)

tan_mu = 0.0
tan_std = 0.9
tan_xbar = exp(tan_mu + tan_std^2 / 2)
tan_xstd = sqrt((exp(tan_std^2) - 1) * exp(2 * tan_mu + tan_std^2))
tan_eq10 = -1.01 + 0.11 * ach / 1.0 + 0.33 * sqrt(roomarea) / roomheight + 0.01 * roomheight / 1.0

# fraction of people who spend some time in the near field
for (iday in 1:ndays) {
  fnear[, iday] <- runif(nindex, min = 0.0, max = 0.2)
}

# number of people who spend some time in the near field
nsus_near = round(fnear * nsus)
nsus_near = apply(nsus_near, MARGIN = c(1,2), FUN = as.integer)

fnear = matrix(0.0, nrow = nindex, ncol = ndays)
for (iday in 1:ndays) {
  for (kind in 1:nindex) {
    if (nsus[kind, iday] > 0) {
      fnear[kind, iday] <- nsus_near[kind, iday] / nsus[kind, iday]
    }
  }
}

for (iday in 1:ndays) {
  for (kind in 1:nindex) {
    if (nsus[kind, iday] > 0) {
      tnear[kind, iday] <- tend[kind, iday] / nsus_near[kind, iday]
    }
  }
}

for (iday in 1:ndays) {
  for (ibin in 1:ibinmax) {
    # Aerosol volume emission rate in each bin (mL/s)
    # nemm_talk and nemm_breathe are for lowbreathe  - scale by actual breathing
    vaerosol <- talkfraction[, iday] * (nemm_talk * vdotbreathe[, iday] / lowbreathe) * vemtalk[ibin] + (1.0 - talkfraction[, iday]) * (nemm_breathe * vdotbreathe[, iday] / lowbreathe) * vembreathe[ibin]
    vaerosoltalk <- nemm_talk * vdotbreathe[, iday] / lowbreathe * vemtalk[ibin]
    vaerosolbreathe <- nemm_breathe * vdotbreathe[, iday] / lowbreathe * vembreathe[ibin]

    # Total aerosol volume emission rate summed over all bins (mL/s)
    vaerosoltot[, iday] <- vaerosoltot[, iday] + vaerosol
    vaerosoltottalk[, iday] <- vaerosoltottalk[, iday] + vaerosoltalk
    vaerosoltotbreathe[, iday] <- vaerosoltotbreathe[, iday] + vaerosolbreathe

    # Loss coefficient (1/s)
    kloss = (kdecay + kdep[ibin] + ach[, iday]) / 3600.0
    klossm = (kdecay + kdep[ibin] + achm[, iday]) / 3600.0

    # CONCENTRATION IN WELL-MIXED CASE
    rnaconc_1[ibin, , iday] <- ((1.0 - maskfac_out) * vinf_1[, iday] * vaerosol) / (kloss * roomvol[, iday]) * (1.0 - (1 - exp(-kloss * tend[, iday])) / (kloss * tend[, iday]))
    rnaconc_2[ibin, , iday] <- ((1.0 - maskfac_out) * vinf_2[, iday] * vaerosol) / (kloss * roomvol[, iday]) * (1.0 - (1 - exp(-kloss * tend[, iday])) / (kloss * tend[, iday]))
    rnaconc_3[ibin, , iday] <- ((1.0 - maskfac_out) * vinf_3[, iday] * vaerosol) / (kloss * roomvol[, iday]) * (1.0 - (1 - exp(-kloss * tend[, iday])) / (kloss * tend[, iday]))
    rnaconc_4[ibin, , iday] <- ((1.0 - maskfac_out) * vinf_4[, iday] * vaerosol) / (kloss * roomvol[, iday]) * (1.0 - (1 - exp(-kloss * tend[, iday])) / (kloss * tend[, iday]))
    rnaconcm_1[ibin, , iday] <- ((1.0 - maskfac_out) * vinf_1[, iday] * vaerosol) / (klossm * roomvol[, iday]) * (1.0 - (1 - exp(-klossm * tend[, iday])) / (klossm * tend[, iday]))
    rnaconcm_2[ibin, , iday] <- ((1.0 - maskfac_out) * vinf_2[, iday] * vaerosol) / (klossm * roomvol[, iday]) * (1.0 - (1 - exp(-klossm * tend[, iday])) / (klossm * tend[, iday]))
    rnaconcm_3[ibin, , iday] <- ((1.0 - maskfac_out) * vinf_3[, iday] * vaerosol) / (klossm * roomvol[, iday]) * (1.0 - (1 - exp(-klossm * tend[, iday])) / (klossm * tend[, iday]))
    rnaconcm_4[ibin, , iday] <- ((1.0 - maskfac_out) * vinf_4[, iday] * vaerosol) / (klossm * roomvol[, iday]) * (1.0 - (1 - exp(-klossm * tend[, iday])) / (klossm * tend[, iday]))

    # FAR-FIELD CONCENTRATION IN NEAR+FAR CASE
    xrnaconc_far_1[ibin, , iday] <- ((1.0 - maskfac_out) * xvinf_1[, iday] * vaerosol) / (kloss * roomvol[, iday]) * (1.0 - (1 - exp(-kloss * tend[, iday])) / (kloss * tend[, iday]))
    xrnaconc_far_2[ibin, , iday] <- ((1.0 - maskfac_out) * xvinf_2[, iday] * vaerosol) / (kloss * roomvol[, iday]) * (1.0 - (1 - exp(-kloss * tend[, iday])) / (kloss * tend[, iday]))
    xrnaconc_far_3[ibin, , iday] <- ((1.0 - maskfac_out) * xvinf_3[, iday] * vaerosol) / (kloss * roomvol[, iday]) * (1.0 - (1 - exp(-kloss * tend[, iday])) / (kloss * tend[, iday]))
    xrnaconc_far_4[ibin, , iday] <- ((1.0 - maskfac_out) * xvinf_4[, iday] * vaerosol) / (kloss * roomvol[, iday]) * (1.0 - (1 - exp(-kloss * tend[, iday])) / (kloss * tend[, iday]))
    xrnaconcm_far_1[ibin, , iday] <- ((1.0 - maskfac_out) * xvinf_1[, iday] * vaerosol) / (klossm * roomvol[, iday]) * (1.0 - (1 - exp(-klossm * tend[, iday])) / (klossm * tend[, iday]))
    xrnaconcm_far_2[ibin, , iday] <- ((1.0 - maskfac_out) * xvinf_2[, iday] * vaerosol) / (klossm * roomvol[, iday]) * (1.0 - (1 - exp(-klossm * tend[, iday])) / (klossm * tend[, iday]))
    xrnaconcm_far_3[ibin, , iday] <- ((1.0 - maskfac_out) * xvinf_3[, iday] * vaerosol) / (klossm * roomvol[, iday]) * (1.0 - (1 - exp(-klossm * tend[, iday])) / (klossm * tend[, iday]))
    xrnaconcm_far_4[ibin, , iday] <- ((1.0 - maskfac_out) * xvinf_4[, iday] * vaerosol) / (klossm * roomvol[, iday]) * (1.0 - (1 - exp(-klossm * tend[, iday])) / (klossm * tend[, iday]))

    # NEAR-FIELD CONCENTRATION IN NEAR+FAR CASE
    xrnaconc_near_1[ibin, , iday] <- (xdil_1 - 1) / xdil_1 * xrnaconc_far_1[ibin, , iday] + 1.0 / xdil_1 * (1 - maskfac_out) * xvinf_1[, iday] * vaerosol / vdotbreathe[, iday]
    xrnaconc_near_2[ibin, , iday] <- (xdil_2 - 1) / xdil_2 * xrnaconc_far_2[ibin, , iday] + 1.0 / xdil_2 * (1 - maskfac_out) * xvinf_2[, iday] * vaerosol / vdotbreathe[, iday]
    xrnaconc_near_3[ibin, , iday] <- (xdil_3 - 1) / xdil_3 * xrnaconc_far_3[ibin, , iday] + 1.0 / xdil_3 * (1 - maskfac_out) * xvinf_3[, iday] * vaerosol / vdotbreathe[, iday]
    xrnaconc_near_4[ibin, , iday] <- (xdil_4 - 1) / xdil_4 * xrnaconc_far_4[ibin, , iday] + 1.0 / xdil_4 * (1 - maskfac_out) * xvinf_4[, iday] * vaerosol / vdotbreathe[, iday]
    xrnaconcm_near_1[ibin, , iday] <- (xdil_1 - 1) / xdil_1 * xrnaconcm_far_1[ibin, , iday] + 1.0 / xdil_1 * (1 - maskfac_out) * xvinf_1[, iday] * vaerosol / vdotbreathe[, iday]
    xrnaconcm_near_2[ibin, , iday] <- (xdil_2 - 1) / xdil_2 * xrnaconcm_far_2[ibin, , iday] + 1.0 / xdil_2 * (1 - maskfac_out) * xvinf_2[, iday] * vaerosol / vdotbreathe[, iday]
    xrnaconcm_near_3[ibin, , iday] <- (xdil_3 - 1) / xdil_3 * xrnaconcm_far_3[ibin, , iday] + 1.0 / xdil_3 * (1 - maskfac_out) * xvinf_3[, iday] * vaerosol / vdotbreathe[, iday]
    xrnaconcm_near_4[ibin, , iday] <- (xdil_4 - 1) / xdil_4 * xrnaconcm_far_4[ibin, , iday] + 1.0 / xdil_4 * (1 - maskfac_out) * xvinf_4[, iday] * vaerosol / vdotbreathe[, iday]
  }

  # CALCULATE DOSE FOR WELL MIXED CASE (virions/m3)

  for (ibin in 1:ibinmax) {

  }
}





