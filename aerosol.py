###########################################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as colors

from numpy.random import default_rng

from scipy.special import gamma, factorial
from scipy.stats import truncnorm
from scipy.stats import norm
from scipy.stats import gmean
import scipy.stats as stats

from scipy.stats import nbinom
from scipy.stats import exponweib
from scipy.optimize import fmin
from scipy.optimize import curve_fit

import statsmodels.api as sm

import lmfit
from sklearn.metrics import r2_score
import pylab

import pandas as pd

############

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 16
plt.rcParams["axes.linewidth"] = 0.5
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.titlesize"] = 16

############

# Set random seed so results are reproducible

rng = default_rng(190319381062290118519841602886838153481)
float_epsilon = np.finfo(float).eps

############
# FIXED PARAMETERS

# Number of index case
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
daynumber = np.arange(ndays)

# Number of hourly events per day
nhours = 1

# Height of room (m)
roomheight = 3.0

# Decay rate constant (1/hr); use placeholder of 1/hour
kdecay = 1.0

###################################################################################
# RESPIRATOR VOLUME EMISSION RATE

# Define dry bins from 0.1 to 100 um
# Integrate from 0.1 micron to 100 microns
nbins = 30
nleft = -1
nright = 2
dlogd = (nright-nleft)/nbins

dfull = np.zeros(nbins)

dedge = np.logspace(nleft, nright, nbins+1)

for i in range(nbins):
 dfull[i] = np.sqrt(dedge[i+1] * dedge[i])

# Cappa 2 mode model for breathing - note cn here is particles/s
# Note that cn here is from Cappa 3/10/2022 4:02 am email
nmodesb  = 2
cnb = np.array([1.33/(1.33+0.88), 0.88/(1.33+0.88)])
cmdb = np.array([0.45, 0.87])
gsdb = np.array([1.58,1.47])

# Cappa 3 mode model for talking - note cn here is particles/s
# Note that cn here is from Cappa 2/27/2022 2:52 pm email
nmodest  = 3
cnt = np.array([4.3/(4.3+0.1+0.011), 0.1/(4.3+0.1+0.011), 0.011/(4.3+0.1+0.011)])
cmdt = np.array([0.45, 3.3, 10.5])
gsdt = np.array([1.96,1.6, 1.6])

# Set max diameter to 25 um dry ie 100 um wet
# ibinmax = 24
# Set max diameter to 5 um dry ie 20 um wet
ibinmax = 17
d = dfull[0:ibinmax]
nemb = np.zeros((ibinmax,nmodesb))
vemb = np.zeros((ibinmax,nmodesb))
nemt = np.zeros((ibinmax,nmodest))
vemt = np.zeros((ibinmax,nmodest))

# Calculate number (1/s) and volume (um3/s) in each size bin for breathing for 1 p/s
for im in range(nmodesb):
 for i in range(ibinmax):
     nemb[i,im] = dlogd * np.log(10.0) * cnb[im]/(np.sqrt(2.0*np.pi)*np.log(gsdb[im])) * np.exp(-((np.log(d[i]) - np.log(cmdb[im]))**2)/(2.0*(np.log(gsdb[im]))**2))
     vemb[i,im] = nemb[i,im] * (np.pi * (d[i])**3)/6.0

# Calculate number (1/s) and volume (um3/s) in each size bin for talking for 1 p/s
for im in range(nmodest):
 for i in range(ibinmax):
     nemt[i,im] = dlogd * np.log(10.0) * cnt[im]/(np.sqrt(2.0*np.pi)*np.log(gsdt[im])) * np.exp(-((np.log(d[i]) - np.log(cmdt[im]))**2)/(2.0*(np.log(gsdt[im]))**2))
     vemt[i,im] = nemt[i,im] * (np.pi * (d[i])**3)/6.0

# Total wet volume emission rate for each bin for breathing for p/s = 1.0 in mL/s
# Scale by 64 for dry to wet volume
# Multiply by 1e-12 to convert um3/s to mL/s
vembreathe = 64.0 * 1.0e-12 * np.sum(vemb,axis=1)

# Total wet volume emission rate for for each bin talking for p/s = 1.0 in mL/s
# Scale by 64 for dry to wet volume
# Multiply by 1e-12 to convert um3/s to mL/s
vemtalk = 64.0 * 1.0e-12 * np.sum(vemt,axis=1)

# Specify distribution of p/s for breathing
# From Cappa's 3/11 and 3/12 emails (also see my email to him on 3/17 at 2.40 pm) for breathing
# Scale by factor of 8 (ie add ln(8) to mu) to account for APS sampling
mub = np.log(8) - 0.48*np.log(10)
sigmab = 0.85*np.log(10)/(np.sqrt(2))
nemm_breathe = rng.lognormal(mub, sigmab, size=nindex)

# Specify distribution of p/s for talking
# From Asadi, fig 4B, draw samples of number of particles/s emitted for talking
# Mean = 3.4 particles/s, Std = 2.7 particles/s
# Scale by factor of 13 (ie multiply mean and std dev by 13) to account for APS sampling
# Truncate log-normal at 260 particles/s (limit of 20 from Fig 5A in Asadi et al. scaled by factor of 13)
mean = 3.4*13
std = 2.7*13
mu = np.log(mean*mean/np.sqrt(mean*mean + std*std))
sigma = np.sqrt(np.log(1.0 + std*std/(mean*mean)))
lmax = np.log(20.0*13)
bclip = (lmax - mu)/sigma
lclip = truncnorm.rvs(-np.inf, bclip, loc=mu, scale=sigma, size=nindex,random_state=rng)
nemm_talk = np.exp(lclip)

# Set p/s from breathing to no more than 10% from talking
for i in range(nindex):
    if nemm_breathe[i] > 0.1*nemm_talk[i] :
       nemm_breathe[i] = 0.1*nemm_talk[i]

# Set p/s from breathing to a max of 24 p/s (limit of 3 from Fig 5A in Asadi et al. scaled by factor of 8)
nemm_breathe = np.clip(nemm_breathe, None, 24)

print('nemm_talk',np.mean(nemm_talk),np.median(nemm_talk),np.std(nemm_talk),np.min(nemm_talk),np.max(nemm_talk))
print('nemm_breathe',np.mean(nemm_breathe),np.median(nemm_breathe),np.std(nemm_breathe),np.min(nemm_breathe),np.max(nemm_breathe))

vbreatheindex = np.zeros(nindex)
vtalkindex = np.zeros(nindex)

# Wet volume emission rate (mL/s) for talking and breathing for each index case
for n in range(nindex):
 vbreatheindex[n] = nemm_breathe[n] * np.sum(vembreathe)
 vtalkindex[n] = nemm_talk[n] * np.sum(vemtalk)

vbreatheindex_median = np.median(vbreatheindex)
vtalkindex_median = np.median(vtalkindex)
vbreatheindex_mean = np.mean(vbreatheindex)
vtalkindex_mean = np.mean(vtalkindex)

vtalkindex_75 = np.percentile(vtalkindex, 75)
vtalkindex_90 = np.percentile(vtalkindex, 90)

vratio = vtalkindex/vbreatheindex
############
# BREATHING RATE

# Breathing rate in m3/s
# Uniform distribution from sitting to light activity from Table 3 Henriques et al. (2022)

vdotbreathe = np.zeros((nindex,ndays))

lowbreathe = 0.51/3600.0
highbreathe = 1.24/3600.0
# lowbreathe = 1.0e-4
# highbreathe = 2.0e-4

# Implicit assumption here is that breathing rate is same for all individuals in room
br_activity = [0.51/3600.0, 0.57/3600.0, 1.24/3600]
activityindex = np.zeros((nindex,ndays), dtype=int)
for iday in range(ndays):
    activityindex[:,iday] = rng.choice(3,nindex)
    for iindex in range(nindex):
        vdotbreathe[iindex,iday]=br_activity[activityindex[iindex,iday]]

############
# VIRAL LOAD AND INFECTIOUSNESS

# Viral load NP (copies)
vnose = np.zeros((nindex,ndays))
kevn = np.loadtxt('data/Ke_Log10VLs_10000.csv', delimiter = ',', skiprows = 1)
# Multiply vnose by 3 to account for 3ml of VMT 
for iday in range(ndays):
    irow = 10*iday+5
    vnose[:,iday] = (10**kevn[irow,1:10001])*3

# Viral load saliva (copies/mL saliva)
vsalv = np.zeros((nindex,ndays))
kevs = np.loadtxt('data/Ke_Log10VLs_saliva_10000.csv', delimiter = ',', skiprows = 1)
for iday in range(ndays):
    irow = 10*iday+5
    vsalv[:,iday] = (10**kevs[irow,1:10001])

# Infectiousness (arbitrary units)
vinfectiousness = np.zeros((nindex,ndays))
kevi = np.loadtxt('data/Ke_Infectiousness_10000.csv', delimiter = ',', skiprows = 1)
# Multiply vnose by 3 to account for 3ml of VMT 
for iday in range(ndays):
    irow = 10*iday+5
    vinfectiousness[:,iday] = 3.0*kevi[irow,1:10001]

##############
# ACH

ach = np.zeros((nindex,ndays))
achm = np.zeros((nindex,ndays))

# From Persily et al. (2005)
meanach = 2.00
stdach = 2.45

mu_ach = np.log(meanach*meanach/np.sqrt(meanach*meanach + stdach*stdach))
sigma_ach = np.sqrt(np.log(1.0 + stdach*stdach/(meanach*meanach)))

for iday in range(ndays):
    ach[:,iday] = rng.lognormal(mu_ach,sigma_ach, size=nindex)
    achm[:, iday] = np.clip(ach[:, iday], achlimit, None)

##############
# DEPOSITION

# Calculate diameter-dependent deposition velocity based on Lai and Nazaroff (2000)
#
# Air viscosity (kg/m-s); Henriques (2022)
mu_air = 1.8e-5

# Air mean free path (um); Equation 8.71, Seinfeld and Pandis (1998)
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

# Slip correction factor; Equation 8.34, Seinfeld and Pandis (1998)
Cc = 1.0 + Kn * (1.257 + 0.4 * np.exp(-1.1/Kn))

# Particle settling velocity (m/s); Equation 8.42, Seinfeld and Pandis (1998)
setvel_particle = ((d * 1.0e-6)**2) * (dens_particle - dens_air) * grav * Cc /(18.0 * mu_air)  

# Friction velocity (m/s); assumed from Lai and Nazaroff (2000)
ustar = 0.01

# Air kinematic viscosity (m2/s)
kinem_air = mu_air/dens_air  
#
# Boltzman constant (m2 kg/s2-K)
boltz = 1.38e-23 

# temperature (K)
temperature = 298.0

# Particle diffusivity in air (m2/s); Equation 8.73, Seinfeld and Pandis (1998)
diff_particle = boltz * temperature * Cc /(3 * np.pi * mu_air * d*1.0e-6)

# Schmidt number
Sc  = kinem_air/diff_particle

# Calculate deposition velocity (m/s) for upward facing horizontal surface; Table 2, Lai and Nazaroff (2000)

# rplus; Lai and Nazaroff (2000)
rplus = d*1.0e-6/2.0 * ustar/kinem_air 

scpow = Sc**(-1.0/3.0)

aconst = 0.5 * np.log(((10.92*scpow + 4.3)**3)/(1.0/Sc + 0.0609)) + np.sqrt(3) * np.arctan((8.6 - 10.92*scpow)/(np.sqrt(3) * 10.92 * scpow))

bconst = 0.5 * np.log(((10.92*scpow + rplus)**3)/(1.0/Sc + 7.669e-4*(rplus**3))) + np.sqrt(3) * np.arctan((2.0*rplus - 10.92*scpow)/(np.sqrt(3) * 10.92 * scpow))

integ = 3.64 * (Sc**0.66667) * (aconst - bconst) + 39.0

# Deposition velocity (m/s)
dep = setvel_particle/(1.0 - np.exp(-setvel_particle*integ/ustar))

# Deposition rate coefficient (1/hr)
kdep = 3600.0*dep/roomheight

##############
# DIAMETER-DEPENDENT DEPOSITION FRACTION IN RESPIRATORY SYSTEM - from Henriques et al. (2022)

ifrac = 1.0 - 0.5*(1.0 - 1.0/(1.0 + 0.00076*(d**2.8)))    
fdep = ifrac * (0.0587 + 0.911/(1.0 + np.exp(4.77 + 1.485*np.log(d))) + 0.943/(1.0 + np.exp(0.508 - 2.58*np.log(d))))

##############

# AREA PER-PERSON (m**2/person)

areapp = np.zeros((nindex,ndays))
densx = np.zeros((nindex,ndays))

# ASHRAE default densities
for iday in range(ndays):
    areapp[:,iday] = rng.uniform(1.0,10.0,size=nindex)
    densx[:,iday] = 100.0/areapp[:,iday]

######
# NUMBER OF SUSCEPTIBLE CONTACTS

nsus = np.zeros((nindex,ndays))

for iday in range(ndays):
    nsus[:,iday] = np.rint(rng.gamma(20.0/30.0,30,size=nindex))

nsus = nsus.astype(np.int_)

###############
# EVENT DURATION (s)

tend = np.zeros((nindex,ndays))

mu_tend = np.log(1.0)
sigma_tend=0.5
lntend_min = np.log(0.1)
lntend_max = np.log(3.0+float_epsilon)
aclip = (lntend_min - mu_tend)/sigma_tend
bclip = (lntend_max - mu_tend)/sigma_tend
for iday in range(ndays):
    lntend_clip = truncnorm.rvs(aclip, bclip, loc=mu_tend, scale=sigma_tend, size=nindex,random_state=rng)
    tend[:, iday] = 3600.0*np.exp(lntend_clip)

###############
# TALK TIME FRACTION

talkfraction = np.zeros((nindex,ndays))

for iday in range(ndays):
    talkfraction[:,iday] = rng.uniform(low=0.05,high=0.25+float_epsilon, size=nindex)

###############
# CALCULATE AREA OF ROOM (m**2)
roomarea = (1 + nsus)*areapp

# CALCULATE VOLUME OF ROOM (m**3)

roomvol = roomheight*roomarea

###############
# CALCULATE VOLUME FLOW RATE (m3/s)

volflow = roomvol*ach/3600.0
volflowmin = 1.e-3 * ( 5.0 * (1 + nsus) + 1.0 * roomarea)
volflownew = np.clip(volflow, volflowmin, None)
achm = volflownew*3600.0/roomvol
volflowpp = volflow/(1+nsus)
volflowppm = volflownew/(1+nsus)

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
hpower_3 = rng.uniform(low=0.5, high=0.75, size = nindex)
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
vinf_1 = const_1*(vsalv**hpower_1)
vinf_2 = const_2*(vnose**hpower_2)
vinf_3 = np.zeros((nindex,ndays))
vinf_4 = np.zeros((nindex,ndays))
for iday in range(ndays):
    vinf_3[:,iday]  = const_3*(vsalv[:,iday]**hpower_3[:])
    vinf_4[:,iday]  = const_4*(vnose[:,iday]**hpower_4[:])

xvinf_1 = xconst_1*(vsalv**hpower_1)
xvinf_2 = xconst_2*(vnose**hpower_2)
xvinf_3 = np.zeros((nindex,ndays))
xvinf_4 = np.zeros((nindex,ndays))
for iday in range(ndays):
    xvinf_3[:,iday]  = xconst_3*(vsalv[:,iday]**hpower_3[:])
    xvinf_4[:,iday]  = xconst_4*(vnose[:,iday]**hpower_4[:])


# ID50 in virions
# SET TO 1 to INCLUDE IN SCALING CONSTANT
id50 = 1.0

# ARRAYS FOR WELL-MIXED CASE
rnaconc_1 = np.zeros((ibinmax,nindex,ndays))
rnaconc_2 = np.zeros((ibinmax,nindex,ndays))
rnaconc_3 = np.zeros((ibinmax,nindex,ndays))
rnaconc_4 = np.zeros((ibinmax,nindex,ndays))
dose_1 = np.zeros((nindex,ndays))
dose_2 = np.zeros((nindex,ndays))
dose_3 = np.zeros((nindex,ndays))
dose_4 = np.zeros((nindex,ndays))
pinf_1 = np.zeros((nindex,ndays))
pinf_2 = np.zeros((nindex,ndays))
pinf_3 = np.zeros((nindex,ndays))
pinf_4 = np.zeros((nindex,ndays))
ninf_1 = np.zeros((nindex,ndays))
ninf_2 = np.zeros((nindex,ndays))
ninf_3 = np.zeros((nindex,ndays))
ninf_4 = np.zeros((nindex,ndays))
ninf_1x = np.zeros((nindex,ndays))
ninf_2x = np.zeros((nindex,ndays))
ninf_3x = np.zeros((nindex,ndays))
ninf_4x = np.zeros((nindex,ndays))
Hr = np.zeros((nindex,ndays))
HrV = np.zeros((nindex,ndays))
Hrtrue = np.zeros((nindex,ndays))
vbreathetrue = np.zeros((nindex,ndays))
vtalktrue = np.zeros((nindex,ndays))

rnaconcm_1 = np.zeros((ibinmax,nindex,ndays))
rnaconcm_2 = np.zeros((ibinmax,nindex,ndays))
rnaconcm_3 = np.zeros((ibinmax,nindex,ndays))
rnaconcm_4 = np.zeros((ibinmax,nindex,ndays))
dosem_1 = np.zeros((nindex,ndays))
dosem_2 = np.zeros((nindex,ndays))
dosem_3 = np.zeros((nindex,ndays))
dosem_4 = np.zeros((nindex,ndays))
pinfm_1 = np.zeros((nindex,ndays))
pinfm_2 = np.zeros((nindex,ndays))
pinfm_3 = np.zeros((nindex,ndays))
pinfm_4 = np.zeros((nindex,ndays))
ninfm_1 = np.zeros((nindex,ndays))
ninfm_2 = np.zeros((nindex,ndays))
ninfm_3 = np.zeros((nindex,ndays))
ninfm_4 = np.zeros((nindex,ndays))
ninfm_1x = np.zeros((nindex,ndays))
ninfm_2x = np.zeros((nindex,ndays))
ninfm_3x = np.zeros((nindex,ndays))
ninfm_4x = np.zeros((nindex,ndays))

# ARRAYS FOR NEAR+FAR CASE

nsus_near = np.zeros((nindex,ndays))
fnear = np.zeros((nindex,ndays))
tnear = np.zeros((nindex,ndays))

xrnaconc_far_1 = np.zeros((ibinmax,nindex,ndays))
xrnaconc_far_2 = np.zeros((ibinmax,nindex,ndays))
xrnaconc_far_3 = np.zeros((ibinmax,nindex,ndays))
xrnaconc_far_4 = np.zeros((ibinmax,nindex,ndays))
xrnaconc_near_1 = np.zeros((ibinmax,nindex,ndays))
xrnaconc_near_2 = np.zeros((ibinmax,nindex,ndays))
xrnaconc_near_3 = np.zeros((ibinmax,nindex,ndays))
xrnaconc_near_4 = np.zeros((ibinmax,nindex,ndays))
xdose_far_1 = np.zeros((nindex,ndays))
xdose_far_2 = np.zeros((nindex,ndays))
xdose_far_3 = np.zeros((nindex,ndays))
xdose_far_4 = np.zeros((nindex,ndays))
xdose_nearfar_1 = np.zeros((nindex,ndays))
xdose_nearfar_2 = np.zeros((nindex,ndays))
xdose_nearfar_3 = np.zeros((nindex,ndays))
xdose_nearfar_4 = np.zeros((nindex,ndays))
xpinf_far_1 = np.zeros((nindex,ndays))
xpinf_far_2 = np.zeros((nindex,ndays))
xpinf_far_3 = np.zeros((nindex,ndays))
xpinf_far_4 = np.zeros((nindex,ndays))
xpinf_nearfar_1 = np.zeros((nindex,ndays))
xpinf_nearfar_2 = np.zeros((nindex,ndays))
xpinf_nearfar_3 = np.zeros((nindex,ndays))
xpinf_nearfar_4 = np.zeros((nindex,ndays))
xpinf_1 = np.zeros((nindex,ndays))
xpinf_2 = np.zeros((nindex,ndays))
xpinf_3 = np.zeros((nindex,ndays))
xpinf_4 = np.zeros((nindex,ndays))
xninf_1 = np.zeros((nindex,ndays))
xninf_2 = np.zeros((nindex,ndays))
xninf_3 = np.zeros((nindex,ndays))
xninf_4 = np.zeros((nindex,ndays))

xrnaconcm_far_1 = np.zeros((ibinmax,nindex,ndays))
xrnaconcm_far_2 = np.zeros((ibinmax,nindex,ndays))
xrnaconcm_far_3 = np.zeros((ibinmax,nindex,ndays))
xrnaconcm_far_4 = np.zeros((ibinmax,nindex,ndays))
xrnaconcm_near_1 = np.zeros((ibinmax,nindex,ndays))
xrnaconcm_near_2 = np.zeros((ibinmax,nindex,ndays))
xrnaconcm_near_3 = np.zeros((ibinmax,nindex,ndays))
xrnaconcm_near_4 = np.zeros((ibinmax,nindex,ndays))
xdosem_far_1 = np.zeros((nindex,ndays))
xdosem_far_2 = np.zeros((nindex,ndays))
xdosem_far_3 = np.zeros((nindex,ndays))
xdosem_far_4 = np.zeros((nindex,ndays))
xdosem_nearfar_1 = np.zeros((nindex,ndays))
xdosem_nearfar_2 = np.zeros((nindex,ndays))
xdosem_nearfar_3 = np.zeros((nindex,ndays))
xdosem_nearfar_4 = np.zeros((nindex,ndays))
xpinfm_far_1 = np.zeros((nindex,ndays))
xpinfm_far_2 = np.zeros((nindex,ndays))
xpinfm_far_3 = np.zeros((nindex,ndays))
xpinfm_far_4 = np.zeros((nindex,ndays))
xpinfm_nearfar_1 = np.zeros((nindex,ndays))
xpinfm_nearfar_2 = np.zeros((nindex,ndays))
xpinfm_nearfar_3 = np.zeros((nindex,ndays))
xpinfm_nearfar_4 = np.zeros((nindex,ndays))
xpinfm_1 = np.zeros((nindex,ndays))
xpinfm_2 = np.zeros((nindex,ndays))
xpinfm_3 = np.zeros((nindex,ndays))
xpinfm_4 = np.zeros((nindex,ndays))
xninfm_1 = np.zeros((nindex,ndays))
xninfm_2 = np.zeros((nindex,ndays))
xninfm_3 = np.zeros((nindex,ndays))
xninfm_4 = np.zeros((nindex,ndays))

Hrm = np.zeros((nindex,ndays))

vaerosoltot = np.zeros((nindex,ndays))
vaerosoltottalk = np.zeros((nindex,ndays))
vaerosoltotbreathe = np.zeros((nindex,ndays))

tan_mu = 0.0
tan_std = 0.9
tan_xbar = np.exp(tan_mu + tan_std*tan_std/2)
tan_xstd = np.sqrt((np.exp(tan_std*tan_std) - 1) * np.exp(2.0*tan_mu + tan_std*tan_std))
tan_eq10 = -1.01 + 0.11*ach/1.0 + 0.33*np.sqrt(roomarea)/roomheight + 0.01*roomheight/1.0


# fraction of people who spend some time in the near field
for iday in range(ndays):
    fnear[:,iday] = rng.uniform(low=0.0, high=0.2, size = nindex)

# number of people who spend some time in the near field
nsus_near = np.rint(fnear * nsus)
nsus_near = nsus_near.astype(np.int_)

fnear[:,:] = 0.0
for iday in range(ndays):
    for kind in range(nindex):
        if nsus[kind,iday] > 0:
           fnear[kind,iday] = nsus_near[kind,iday]/nsus[kind,iday]

for iday in range(ndays):
    for kind in range(nindex):
        if nsus_near[kind,iday] > 0:
           tnear[kind,iday] = tend[kind,iday]/nsus_near[kind,iday]

for iday in range(ndays):
    for ibin in range(ibinmax):

        # Aerosol volume emission rate in each bin (mL/s)
        # nemm_talk and nemm_breathe are for lowbreathe  - scale by actual breathing
        vaerosol = talkfraction[:,iday] * (nemm_talk[:] * vdotbreathe[:,iday]/lowbreathe) * vemtalk[ibin] + (1.0 - talkfraction[:,iday]) * (nemm_breathe[:] * vdotbreathe[:,iday]/lowbreathe) * vembreathe[ibin] 
        vaerosoltalk = nemm_talk[:] * vdotbreathe[:,iday]/lowbreathe * vemtalk[ibin] 
        vaerosolbreathe = nemm_breathe[:] * vdotbreathe[:,iday]/lowbreathe * vembreathe[ibin]

        # Total aerosol volume emission rate summed over all bins (mL/s)
        vaerosoltot[:,iday] = vaerosoltot[:,iday] + vaerosol[:]
        vaerosoltottalk[:,iday] = vaerosoltottalk[:,iday] + vaerosoltalk[:]
        vaerosoltotbreathe[:,iday] = vaerosoltotbreathe[:,iday] + vaerosolbreathe[:]

        # Loss coefficient (1/s)
        kloss = (kdecay + kdep[ibin] + ach[:,iday])/3600.0 
        klossm = (kdecay + kdep[ibin] + achm[:,iday])/3600.0 

        # CONCENTRATION IN WELL-MIXED CASE
        rnaconc_1[ibin,:,iday] = ((1.0 - maskfac_out)*vinf_1[:,iday]*vaerosol[:])/(kloss[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-kloss[:] * tend[:,iday]))/(kloss[:] * tend[:,iday]))
        rnaconc_2[ibin,:,iday] = ((1.0 - maskfac_out)*vinf_2[:,iday]*vaerosol[:])/(kloss[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-kloss[:] * tend[:,iday]))/(kloss[:] * tend[:,iday]))
        rnaconc_3[ibin,:,iday] = ((1.0 - maskfac_out)*vinf_3[:,iday]*vaerosol[:])/(kloss[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-kloss[:] * tend[:,iday]))/(kloss[:] * tend[:,iday]))
        rnaconc_4[ibin,:,iday] = ((1.0 - maskfac_out)*vinf_4[:,iday]*vaerosol[:])/(kloss[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-kloss[:] * tend[:,iday]))/(kloss[:] * tend[:,iday]))
        rnaconcm_1[ibin,:,iday] = ((1.0 - maskfac_out)*vinf_1[:,iday]*vaerosol[:])/(klossm[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-klossm[:] * tend[:,iday]))/(klossm[:] * tend[:,iday]))
        rnaconcm_2[ibin,:,iday] = ((1.0 - maskfac_out)*vinf_2[:,iday]*vaerosol[:])/(klossm[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-klossm[:] * tend[:,iday]))/(klossm[:] * tend[:,iday]))
        rnaconcm_3[ibin,:,iday] = ((1.0 - maskfac_out)*vinf_3[:,iday]*vaerosol[:])/(klossm[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-klossm[:] * tend[:,iday]))/(klossm[:] * tend[:,iday]))
        rnaconcm_4[ibin,:,iday] = ((1.0 - maskfac_out)*vinf_4[:,iday]*vaerosol[:])/(klossm[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-klossm[:] * tend[:,iday]))/(klossm[:] * tend[:,iday]))
        
        # FAR-FIELD CONCENTRATION IN NEAR+FAR CASE
        xrnaconc_far_1[ibin,:,iday] = ((1.0 - maskfac_out)*xvinf_1[:,iday]*vaerosol[:])/(kloss[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-kloss[:] * tend[:,iday]))/(kloss[:] * tend[:,iday]))
        xrnaconc_far_2[ibin,:,iday] = ((1.0 - maskfac_out)*xvinf_2[:,iday]*vaerosol[:])/(kloss[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-kloss[:] * tend[:,iday]))/(kloss[:] * tend[:,iday]))
        xrnaconc_far_3[ibin,:,iday] = ((1.0 - maskfac_out)*xvinf_3[:,iday]*vaerosol[:])/(kloss[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-kloss[:] * tend[:,iday]))/(kloss[:] * tend[:,iday]))
        xrnaconc_far_4[ibin,:,iday] = ((1.0 - maskfac_out)*xvinf_4[:,iday]*vaerosol[:])/(kloss[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-kloss[:] * tend[:,iday]))/(kloss[:] * tend[:,iday]))
        xrnaconcm_far_1[ibin,:,iday] = ((1.0 - maskfac_out)*xvinf_1[:,iday]*vaerosol[:])/(klossm[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-klossm[:] * tend[:,iday]))/(klossm[:] * tend[:,iday]))
        xrnaconcm_far_2[ibin,:,iday] = ((1.0 - maskfac_out)*xvinf_2[:,iday]*vaerosol[:])/(klossm[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-klossm[:] * tend[:,iday]))/(klossm[:] * tend[:,iday]))
        xrnaconcm_far_3[ibin,:,iday] = ((1.0 - maskfac_out)*xvinf_3[:,iday]*vaerosol[:])/(klossm[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-klossm[:] * tend[:,iday]))/(klossm[:] * tend[:,iday]))
        xrnaconcm_far_4[ibin,:,iday] = ((1.0 - maskfac_out)*xvinf_4[:,iday]*vaerosol[:])/(klossm[:]*roomvol[:,iday]) * (1.0 - (1 - np.exp(-klossm[:] * tend[:,iday]))/(klossm[:] * tend[:,iday]))

        # NEAR-FIELD CONCENTRATION IN NEAR+FAR CASE
        xrnaconc_near_1[ibin,:,iday] = (xdil_1-1)/xdil_1*xrnaconc_far_1[ibin,:,iday] + 1.0/xdil_1*(1-maskfac_out)*xvinf_1[:,iday]*vaerosol[:]/vdotbreathe[:,iday]
        xrnaconc_near_2[ibin,:,iday] = (xdil_2-1)/xdil_2*xrnaconc_far_2[ibin,:,iday] + 1.0/xdil_2*(1-maskfac_out)*xvinf_2[:,iday]*vaerosol[:]/vdotbreathe[:,iday]
        xrnaconc_near_3[ibin,:,iday] = (xdil_3-1)/xdil_3*xrnaconc_far_3[ibin,:,iday] + 1.0/xdil_3*(1-maskfac_out)*xvinf_3[:,iday]*vaerosol[:]/vdotbreathe[:,iday]
        xrnaconc_near_4[ibin,:,iday] = (xdil_4-1)/xdil_4*xrnaconc_far_4[ibin,:,iday] + 1.0/xdil_4*(1-maskfac_out)*xvinf_4[:,iday]*vaerosol[:]/vdotbreathe[:,iday]
        xrnaconcm_near_1[ibin,:,iday] = (xdil_1-1)/xdil_1*xrnaconcm_far_1[ibin,:,iday] + 1.0/xdil_1*(1-maskfac_out)*xvinf_1[:,iday]*vaerosol[:]/vdotbreathe[:,iday]
        xrnaconcm_near_2[ibin,:,iday] = (xdil_2-1)/xdil_2*xrnaconcm_far_2[ibin,:,iday] + 1.0/xdil_2*(1-maskfac_out)*xvinf_2[:,iday]*vaerosol[:]/vdotbreathe[:,iday]
        xrnaconcm_near_3[ibin,:,iday] = (xdil_3-1)/xdil_3*xrnaconcm_far_3[ibin,:,iday] + 1.0/xdil_3*(1-maskfac_out)*xvinf_3[:,iday]*vaerosol[:]/vdotbreathe[:,iday]
        xrnaconcm_near_4[ibin,:,iday] = (xdil_4-1)/xdil_4*xrnaconcm_far_4[ibin,:,iday] + 1.0/xdil_4*(1-maskfac_out)*xvinf_4[:,iday]*vaerosol[:]/vdotbreathe[:,iday]

    # CALCULATE DOSE FOR WELL MIXED CASE (virions/m3)

    for ibin in range(ibinmax):
        # DOSE FOR SUSCEPTIBLE IN WELL-MIXED CASE
        dose_1[:,iday] = dose_1[:,iday] + (1-maskfac_in)*fdep[ibin]*rnaconc_1[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 
        dose_2[:,iday] = dose_2[:,iday] + (1-maskfac_in)*fdep[ibin]*rnaconc_2[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 
        dose_3[:,iday] = dose_3[:,iday] + (1-maskfac_in)*fdep[ibin]*rnaconc_3[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 
        dose_4[:,iday] = dose_4[:,iday] + (1-maskfac_in)*fdep[ibin]*rnaconc_4[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 
        dosem_1[:,iday] = dosem_1[:,iday] + (1-maskfac_in)*fdep[ibin]*rnaconcm_1[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 
        dosem_2[:,iday] = dosem_2[:,iday] + (1-maskfac_in)*fdep[ibin]*rnaconcm_2[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 
        dosem_3[:,iday] = dosem_3[:,iday] + (1-maskfac_in)*fdep[ibin]*rnaconcm_3[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 
        dosem_4[:,iday] = dosem_4[:,iday] + (1-maskfac_in)*fdep[ibin]*rnaconcm_4[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 

        # DOSE FOR SUSCEPTIBLE WHO ARE ONLY IN FAR-FIELD IN NEAR+FAR CASE
        xdose_far_1[:,iday] = xdose_far_1[:,iday] + (1-maskfac_in)*fdep[ibin]*xrnaconc_far_1[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 
        xdose_far_2[:,iday] = xdose_far_2[:,iday] + (1-maskfac_in)*fdep[ibin]*xrnaconc_far_2[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 
        xdose_far_3[:,iday] = xdose_far_3[:,iday] + (1-maskfac_in)*fdep[ibin]*xrnaconc_far_3[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 
        xdose_far_4[:,iday] = xdose_far_4[:,iday] + (1-maskfac_in)*fdep[ibin]*xrnaconc_far_4[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 
        xdosem_far_1[:,iday] = xdosem_far_1[:,iday] + (1-maskfac_in)*fdep[ibin]*xrnaconcm_far_1[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 
        xdosem_far_2[:,iday] = xdosem_far_2[:,iday] + (1-maskfac_in)*fdep[ibin]*xrnaconcm_far_2[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 
        xdosem_far_3[:,iday] = xdosem_far_3[:,iday] + (1-maskfac_in)*fdep[ibin]*xrnaconcm_far_3[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 
        xdosem_far_4[:,iday] = xdosem_far_4[:,iday] + (1-maskfac_in)*fdep[ibin]*xrnaconcm_far_4[ibin,:,iday] * vdotbreathe[:,iday] * tend[:,iday] 


        # DOSE FOR SUSCEPTIBLE WHO ARE IN FAR FIELD AND NEAR FIELD IN NEAR+FAR CASE
        xdose_nearfar_1[:,iday] = xdose_nearfar_1[:,iday] + (1-maskfac_in)*fdep[ibin]*vdotbreathe[:,iday]*((tend[:,iday]-tnear[:,iday])*xrnaconc_far_1[ibin,:,iday] + xrnaconc_near_1[ibin,:,iday] * tnear[:,iday])  
        xdose_nearfar_2[:,iday] = xdose_nearfar_2[:,iday] + (1-maskfac_in)*fdep[ibin]*vdotbreathe[:,iday]*((tend[:,iday]-tnear[:,iday])*xrnaconc_far_2[ibin,:,iday] + xrnaconc_near_2[ibin,:,iday] * tnear[:,iday])  
        xdose_nearfar_3[:,iday] = xdose_nearfar_3[:,iday] + (1-maskfac_in)*fdep[ibin]*vdotbreathe[:,iday]*((tend[:,iday]-tnear[:,iday])*xrnaconc_far_3[ibin,:,iday] + xrnaconc_near_3[ibin,:,iday] * tnear[:,iday])  
        xdose_nearfar_4[:,iday] = xdose_nearfar_4[:,iday] + (1-maskfac_in)*fdep[ibin]*vdotbreathe[:,iday]*((tend[:,iday]-tnear[:,iday])*xrnaconc_far_4[ibin,:,iday] + xrnaconc_near_4[ibin,:,iday] * tnear[:,iday])  
        xdosem_nearfar_1[:,iday] = xdosem_nearfar_1[:,iday] + (1-maskfac_in)*fdep[ibin]*vdotbreathe[:,iday]*((tend[:,iday]-tnear[:,iday])*xrnaconcm_far_1[ibin,:,iday] + xrnaconcm_near_1[ibin,:,iday] * tnear[:,iday])  
        xdosem_nearfar_2[:,iday] = xdosem_nearfar_2[:,iday] + (1-maskfac_in)*fdep[ibin]*vdotbreathe[:,iday]*((tend[:,iday]-tnear[:,iday])*xrnaconcm_far_2[ibin,:,iday] + xrnaconcm_near_2[ibin,:,iday] * tnear[:,iday])  
        xdosem_nearfar_3[:,iday] = xdosem_nearfar_3[:,iday] + (1-maskfac_in)*fdep[ibin]*vdotbreathe[:,iday]*((tend[:,iday]-tnear[:,iday])*xrnaconcm_far_3[ibin,:,iday] + xrnaconcm_near_3[ibin,:,iday] * tnear[:,iday])  
        xdosem_nearfar_4[:,iday] = xdosem_nearfar_4[:,iday] + (1-maskfac_in)*fdep[ibin]*vdotbreathe[:,iday]*((tend[:,iday]-tnear[:,iday])*xrnaconcm_far_4[ibin,:,iday] + xrnaconcm_near_4[ibin,:,iday] * tnear[:,iday])  

    # CALCULATE INFECTION PROBABILITY
    # NOTE: ID50 IS NOT REQUIRED IN CALCULATION OF *_4 INFECTION PROBABILITIES BECAUSE CORRESPONDING DOSES ARE ALREADY IN SCALED INFECTIOUS UNITS

    # INFECTION PROBABILITY FOR SUSCEPTIBLE IN WELL MIXED CASE
    pinf_1[:,iday] = 1.0 - np.exp(-dose_1[:,iday]/id50)
    pinf_2[:,iday] = 1.0 - np.exp(-dose_2[:,iday]/id50)
    pinf_3[:,iday] = 1.0 - np.exp(-dose_3[:,iday]/id50)
    pinf_4[:,iday] = 1.0 - np.exp(-dose_4[:,iday]/id50)
    pinfm_1[:,iday] = 1.0 - np.exp(-dosem_1[:,iday]/id50)
    pinfm_2[:,iday] = 1.0 - np.exp(-dosem_2[:,iday]/id50)
    pinfm_3[:,iday] = 1.0 - np.exp(-dosem_3[:,iday]/id50)
    pinfm_4[:,iday] = 1.0 - np.exp(-dosem_4[:,iday]/id50)

    # INFECTION PROBABILITY FOR SUSCEPTIBLE WHO ARE ONLY IN FAR FIELD IN THE NEAR+FAR CASE
    xpinf_far_1[:,iday] = 1.0 - np.exp(-xdose_far_1[:,iday]/id50)
    xpinf_far_2[:,iday] = 1.0 - np.exp(-xdose_far_2[:,iday]/id50)
    xpinf_far_3[:,iday] = 1.0 - np.exp(-xdose_far_3[:,iday]/id50)
    xpinf_far_4[:,iday] = 1.0 - np.exp(-xdose_far_4[:,iday]/id50)
    xpinfm_far_1[:,iday] = 1.0 - np.exp(-xdosem_far_1[:,iday]/id50)
    xpinfm_far_2[:,iday] = 1.0 - np.exp(-xdosem_far_2[:,iday]/id50)
    xpinfm_far_3[:,iday] = 1.0 - np.exp(-xdosem_far_3[:,iday]/id50)
    xpinfm_far_4[:,iday] = 1.0 - np.exp(-xdosem_far_4[:,iday]/id50)

    # INFECTION PROBABILITY FOR SUSCEPTIBLE WHO ARE IN FAR FIELD AND NEAR FIELD IN THE NEAR+FAR CASE
    xpinf_nearfar_1[:,iday] = 1.0 - np.exp(-xdose_nearfar_1[:,iday]/id50)
    xpinf_nearfar_2[:,iday] = 1.0 - np.exp(-xdose_nearfar_2[:,iday]/id50)
    xpinf_nearfar_3[:,iday] = 1.0 - np.exp(-xdose_nearfar_3[:,iday]/id50)
    xpinf_nearfar_4[:,iday] = 1.0 - np.exp(-xdose_nearfar_4[:,iday]/id50)
    xpinfm_nearfar_1[:,iday] = 1.0 - np.exp(-xdosem_nearfar_1[:,iday]/id50)
    xpinfm_nearfar_2[:,iday] = 1.0 - np.exp(-xdosem_nearfar_2[:,iday]/id50)
    xpinfm_nearfar_3[:,iday] = 1.0 - np.exp(-xdosem_nearfar_3[:,iday]/id50)
    xpinfm_nearfar_4[:,iday] = 1.0 - np.exp(-xdosem_nearfar_4[:,iday]/id50)

    # INFECTION PROBABILITY FOR NEAR+FAR CASE
    xpinf_1[:,iday] = (1.0-fnear[:,iday])*xpinf_far_1[:,iday] + fnear[:,iday]*xpinf_nearfar_1[:,iday]  
    xpinf_2[:,iday] = (1.0-fnear[:,iday])*xpinf_far_2[:,iday] + fnear[:,iday]*xpinf_nearfar_2[:,iday]  
    xpinf_3[:,iday] = (1.0-fnear[:,iday])*xpinf_far_3[:,iday] + fnear[:,iday]*xpinf_nearfar_3[:,iday]  
    xpinf_4[:,iday] = (1.0-fnear[:,iday])*xpinf_far_4[:,iday] + fnear[:,iday]*xpinf_nearfar_4[:,iday]  
    xpinfm_1[:,iday] = (1.0-fnear[:,iday])*xpinfm_far_1[:,iday] + fnear[:,iday]*xpinfm_nearfar_1[:,iday]  
    xpinfm_2[:,iday] = (1.0-fnear[:,iday])*xpinfm_far_2[:,iday] + fnear[:,iday]*xpinfm_nearfar_2[:,iday]  
    xpinfm_3[:,iday] = (1.0-fnear[:,iday])*xpinfm_far_3[:,iday] + fnear[:,iday]*xpinfm_nearfar_3[:,iday]  
    xpinfm_4[:,iday] = (1.0-fnear[:,iday])*xpinfm_far_4[:,iday] + fnear[:,iday]*xpinfm_nearfar_4[:,iday]  


    # CALCULATE NUMBER OF SECONDARY INFECTIONS

    # NUMBER OF SECONDARY INFECTIONS FOR WELL MIXED CASE
    ninf_1[:,iday] = np.rint(nsus[:,iday]*pinf_1[:,iday])
    ninf_2[:,iday] = np.rint(nsus[:,iday]*pinf_2[:,iday])
    ninf_3[:,iday] = np.rint(nsus[:,iday]*pinf_3[:,iday])
    ninf_4[:,iday] = np.rint(nsus[:,iday]*pinf_4[:,iday])
    ninfm_1[:,iday] = np.rint(nsus[:,iday]*pinfm_1[:,iday])
    ninfm_2[:,iday] = np.rint(nsus[:,iday]*pinfm_2[:,iday])
    ninfm_3[:,iday] = np.rint(nsus[:,iday]*pinfm_3[:,iday])
    ninfm_4[:,iday] = np.rint(nsus[:,iday]*pinfm_4[:,iday])

    # NUBER OF SECONDARY INFECTIONS FOR NEAR+FAR CASE
    xninf_1[:,iday] = np.rint(nsus[:,iday]*xpinf_1[:,iday])
    xninf_2[:,iday] = np.rint(nsus[:,iday]*xpinf_2[:,iday])
    xninf_3[:,iday] = np.rint(nsus[:,iday]*xpinf_3[:,iday])
    xninf_4[:,iday] = np.rint(nsus[:,iday]*xpinf_4[:,iday])
    xninfm_1[:,iday] = np.rint(nsus[:,iday]*xpinfm_1[:,iday])
    xninfm_2[:,iday] = np.rint(nsus[:,iday]*xpinfm_2[:,iday])
    xninfm_3[:,iday] = np.rint(nsus[:,iday]*xpinfm_3[:,iday])
    xninfm_4[:,iday] = np.rint(nsus[:,iday]*xpinfm_4[:,iday])

# CALCULATE R*
print('Number of new infections - WELL-MIXED; h=1.0 SALIVA: ', np.sum(ninf_1))
print('R*: ', np.sum(ninf_1)/nindex)
print('Number of new infections - WELL-MIXED; h=1.0 NOSE: ', np.sum(ninf_2))
print('R*: ', np.sum(ninf_2)/nindex)
print('Number of new infections - WELL-MIXED; h=0.5-0.75 SALIVA: ', np.sum(ninf_3))
print('R*: ', np.sum(ninf_3)/nindex)
print('Number of new infections - WELL-MIXED; h=0.5-0.75 NOSE: ', np.sum(ninf_4))
print('R*: ', np.sum(ninf_4)/nindex)
print('')

print('Number of new infections - NF' + str(xdil_1) + '; h=1.0 SALIVA: ', np.sum(xninf_1))
print('R*: ', np.sum(xninf_1)/nindex)
print('Number of new infections - NF' + str(xdil_2) + '; h=1.0 NOSE: ', np.sum(xninf_2))
print('R*: ', np.sum(xninf_2)/nindex)
print('Number of new infections - NF' + str(xdil_3) + '; h=0.5-0.75 SALIVA: ', np.sum(xninf_3))
print('R*: ', np.sum(xninf_3)/nindex)
print('Number of new infections - NF' + str(xdil_4) + '; h=0.5-0.75 NOSE: ', np.sum(xninf_4))
print('R*: ', np.sum(xninf_4)/nindex)
print('')

print('Number of new infections - SPATIAL VARIABILITY CASE; h=1.0 SALIVA: ', np.sum(ninf_1x))
print('R*: ', np.sum(ninf_1x)/nindex) 
print('Number of new infections - SPATIAL VARIABILITY CASE; h=1.0 NOSE: ', np.sum(ninf_2x))
print('R*: ', np.sum(ninf_2x)/nindex) 
print('Number of new infections - SPATIAL VARIABILITY CASE; h=0.5-0.75 SALIVA: ', np.sum(ninf_3x))
print('R*: ', np.sum(ninf_3x)/nindex) 
print('Number of new infections - SPATIAL VARIABILITY CASE; h=0.5-0.75 NOSE: ', np.sum(ninf_4x))
print('R*: ', np.sum(ninf_4x)/nindex) 
print('')

print('MITIGATION Number of new infections - WELL-MIXED; h=1.0 SALIVA: ', np.sum(ninfm_1))
print('R*: ', np.sum(ninfm_1)/nindex)
print('MITIGATION Number of new infections - WELL-MIXED; h=1.0 NOSE: ', np.sum(ninfm_2))
print('R*: ', np.sum(ninfm_2)/nindex)
print('MITIGATION Number of new infections - WELL-MIXED; h=0.5-0.75 SALIVA: ', np.sum(ninfm_3))
print('R*: ', np.sum(ninfm_3)/nindex)
print('MITIGATION Number of new infections - WELL-MIXED; h=0.5-0.75 NOSE: ', np.sum(ninfm_4))
print('R*: ', np.sum(ninfm_4)/nindex)
print('')

print('MITIGATION Number of new infections - NF' + str(xdil_1) + '; h=1.0 SALIVA: ', np.sum(xninfm_1))
print('MITIGATION R*: ', np.sum(xninfm_1)/nindex) 
print('MITIGATION Number of new infections - NF' + str(xdil_2) + '; h=1.0 NOSE: ', np.sum(xninfm_2))
print('MITIGATION R*: ', np.sum(xninfm_2)/nindex) 
print('MITIGATION Number of new infections - NF' + str(xdil_3) + '; h=0.5-0.75 SALIVA: ', np.sum(xninfm_3))
print('MITIGATION R*: ', np.sum(xninfm_3)/nindex) 
print('MITIGATION Number of new infections - NF' + str(xdil_4) + '; h=0.5-0.75 NOSE: ', np.sum(xninfm_4))
print('MITIGATION R*: ', np.sum(xninfm_4)/nindex) 
print('')

print('MITIGATION Number of new infections - SPATIAL VARIABILITY CASE; h=1.0 SALIVA: ', np.sum(ninfm_1x))
print('MITIGATION R*: ', np.sum(ninfm_1x)/nindex) 
print('MITIGATION Number of new infections - SPATIAL VARIABILITY CASE; h=1.0 NOSE: ', np.sum(ninfm_2x))
print('MITIGATION R*: ', np.sum(ninfm_2x)/nindex) 
print('MITIGATION Number of new infections - SPATIAL VARIABILITY CASE; h=0.5-0.75 SALIVA: ', np.sum(ninfm_3x))
print('MITIGATION R*: ', np.sum(ninfm_3x)/nindex) 
print('MITIGATION Number of new infections - SPATIAL VARIABILITY CASE; h=0.5-0.75 NOSE: ', np.sum(ninfm_4x))
print('MITIGATION R*: ', np.sum(ninfm_4x)/nindex) 
print('')

ninfindex_1 = np.sum(ninf_1, axis=1)
ninfindex_2 = np.sum(ninf_2, axis=1)
ninfindex_3 = np.sum(ninf_3, axis=1)
ninfindex_4 = np.sum(ninf_4, axis=1)
xninfindex_1 = np.sum(xninf_1, axis=1)
xninfindex_2 = np.sum(xninf_2, axis=1)
xninfindex_3 = np.sum(xninf_3, axis=1)
xninfindex_4 = np.sum(xninf_4, axis=1)

nzeros_1 = ninfindex_1[ninfindex_1 == 0]
nzeros_2 = ninfindex_2[ninfindex_2 == 0]
nzeros_3 = ninfindex_3[ninfindex_3 == 0]
nzeros_4 = ninfindex_4[ninfindex_4 == 0]
xnzeros_1 = ninfindex_1[xninfindex_1 == 0]
xnzeros_2 = ninfindex_2[xninfindex_2 == 0]
xnzeros_3 = ninfindex_3[xninfindex_3 == 0]
xnzeros_4 = ninfindex_3[xninfindex_4 == 0]

ninfday_1 = np.sum(ninf_1, axis=0)
ninfday_2 = np.sum(ninf_2, axis=0)
ninfday_3 = np.sum(ninf_3, axis=0)
ninfday_4 = np.sum(ninf_4, axis=0)
xninfday_1 = np.sum(xninf_1, axis=0)
xninfday_2 = np.sum(xninf_2, axis=0)
xninfday_3 = np.sum(xninf_3, axis=0)
xninfday_4 = np.sum(xninf_4, axis=0)

ezeros_1 = ninf_1.flatten()[ninf_1.flatten() == 0]
ezeros_2 = ninf_2.flatten()[ninf_2.flatten() == 0]
ezeros_3 = ninf_3.flatten()[ninf_3.flatten() == 0]
ezeros_4 = ninf_4.flatten()[ninf_4.flatten() == 0]
xezeros_1 = xninf_1.flatten()[xninf_1.flatten() == 0]
xezeros_2 = xninf_2.flatten()[xninf_2.flatten() == 0]
xezeros_3 = xninf_3.flatten()[xninf_3.flatten() == 0]
xezeros_4 = xninf_4.flatten()[xninf_4.flatten() == 0]

print('Max infections/index case WELL-MIXED; h=1.0 SALIVA: ',np.max(ninfindex_1))
print('Max infections/index case WELL-MIXED; h=1.0 NOSE: ',np.max(ninfindex_2))
print('Max infections/index case WELL-MIXED; h=0.5-0.75 SALIVA: ',np.max(ninfindex_3))
print('Max infections/index case WELL-MIXED; h=0.5-0.75 NOSE: ',np.max(ninfindex_4))
print('')

print('Max infections/index case NF' + str(xdil_1) + '; h=1.0 SALIVA: ',np.max(xninfindex_1))
print('Max infections/index case NF' + str(xdil_2) + '; h=1.0 NOSE: ',np.max(xninfindex_2))
print('Max infections/index case NF' + str(xdil_3) + '; h=0.5-0.75 SALIVA: ',np.max(xninfindex_3))
print('Max infections/index case NF' + str(xdil_4) + '; h=0.5-0.75 NOSE: ',np.max(xninfindex_4))
print('')

print('Max infections/event WELL-MIXED; h=1.0 SALIVA: ',np.max(ninf_1))
print('Max infections/event WELL-MIXED; h=1.0 NOSE: ',np.max(ninf_2))
print('Max infections/event WELL-MIXED; h=0.5-0.75 SALIVA: ',np.max(ninf_3))
print('Max infections/event WELL-MIXED; h=0.5-0.75 NOSE: ',np.max(ninf_4))
print('')

print('Max infections/event NF' + str(xdil_1) + '; h=1.0 SALIVA: ',np.max(xninf_1))
print('Max infections/event NF' + str(xdil_2) + '; h=1.0 NOSE: ',np.max(xninf_2))
print('Max infections/event NF' + str(xdil_3) + '; h=0.5-0.75 SALIVA: ',np.max(xninf_3))
print('Max infections/event NF' + str(xdil_4) + '; h=0.5-0.75 NOSE: ',np.max(xninf_4))
print('')

print('Number of index cases with 0 infections WELL-MIXED; h=1.0 SALIVA: ', nzeros_1.size)
print('Number of index cases with 0 infections WELL-MIXED; h=1.0 NOSE: ', nzeros_2.size)
print('Number of index cases with 0 infections WELL-MIXED; h=0.5-0.75 SALIVA: ', nzeros_3.size)
print('Number of index cases with 0 infections WELL-MIXED; h=0.5-0.75 NOSE: ', nzeros_4.size)
print('')

print('Number of index cases with 0 infections NF' + str(xdil_1) + '; h=1.0 SALIVA: ', xnzeros_1.size)
print('Number of index cases with 0 infections NF' + str(xdil_2) + '; h=1.0 NOSE: ', xnzeros_2.size)
print('Number of index cases with 0 infections NF' + str(xdil_3) + '; h=0.5-0.75 SALIVA: ', xnzeros_3.size)
print('Number of index cases with 0 infections NF' + str(xdil_4) + '; h=0.5-0.75 NOSE: ', xnzeros_4.size)
print('')

dzeros_1 = np.zeros(ndays)
dzeros_2 = np.zeros(ndays)
dzeros_3 = np.zeros(ndays)
dzeros_4 = np.zeros(ndays)
xdzeros_1 = np.zeros(ndays)
xdzeros_2 = np.zeros(ndays)
xdzeros_3 = np.zeros(ndays)
xdzeros_4 = np.zeros(ndays)
print('Number of events with 0 infections WELL-MIXED; h=1.0 SALIVA: ', ezeros_1.size)
for iday in range(ndays):
    dzeros_1[iday] = ninf_1[:,iday][ninf_1[:,iday] == 0].size 
    print(iday,ninfday_1[iday],dzeros_1[iday])
print('Number of events with 0 infections WELL-MIXED; h=1.0 NOSE: ', ezeros_2.size)
for iday in range(ndays):
    dzeros_2[iday] = ninf_2[:,iday][ninf_2[:,iday] == 0].size 
    print(iday,ninfday_2[iday],dzeros_2[iday])
print('Number of events with 0 infections WELL-MIXED; h=0.5-0.75 SALIVA: ', ezeros_3.size)
for iday in range(ndays):
    dzeros_3[iday] = ninf_3[:,iday][ninf_3[:,iday] == 0].size 
    print(iday,ninfday_3[iday],dzeros_3[iday])
print('Number of events with 0 infections WELL-MIXED; h=0.5-0.75 NOSE: ', ezeros_4.size)
for iday in range(ndays):
    dzeros_4[iday] = ninf_4[:,iday][ninf_4[:,iday] == 0].size 
    print(iday,ninfday_4[iday],dzeros_4[iday])
print('')

print('Number of events with 0 infections NF' + str(xdil_1) + '; h=1.0 SALIVA: ', xezeros_1.size)
for iday in range(ndays):
    xdzeros_1[iday] = xninf_1[:,iday][xninf_1[:,iday] == 0].size
    print(iday,xninfday_1[iday],xdzeros_1[iday])
print('Number of events with 0 infections NF' + str(xdil_2) + '; h=1.0 NOSE: ', xezeros_2.size)
for iday in range(ndays):
    xdzeros_2[iday] = xninf_2[:,iday][xninf_2[:,iday] == 0].size
    print(iday,xninfday_2[iday],xdzeros_2[iday])
print('Number of events with 0 infections NF' + str(xdil_3) + '; h=0.5-0.75 SALIVA: ', xezeros_3.size)
for iday in range(ndays):
    xdzeros_3[iday] = xninf_3[:,iday][xninf_3[:,iday] == 0].size
    print(iday,xninfday_3[iday],xdzeros_3[iday])
print('Number of events with 0 infections NF' + str(xdil_4) + '; h=0.5-0.75 NOSE: ', xezeros_4.size)
for iday in range(ndays):
    xdzeros_4[iday] = xninf_4[:,iday][xninf_4[:,iday] == 0].size
    print(iday,xninfday_4[iday],xdzeros_4[iday])
print('')
