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
volflowpp = volflow/(1+nsus)

###############
# MAIN CALCULATIONS
# const_1 chosen by trial and error to get R0 ~ 3

# Linear - SALIVA  - ONLY FAR FIELD INTERACTIONS
hpower_1 = 1.0
const_1 = 0.4498

# Infection virions per mL of aerosol volume
vinf_1 = const_1*(vsalv**hpower_1)

# ID50 in virions
# SET TO 1 to INCLUDE IN SCALING CONSTANT
id50 = 1.0

# ARRAYS FOR WELL-MIXED CASE
pinf_1 = np.zeros((nindex,ndays))
ninf_1 = np.zeros((nindex,ndays))
Hr = np.zeros((nindex,ndays))
HrV = np.zeros((nindex,ndays))
Hrtrue = np.zeros((nindex,ndays))
vbreathetrue = np.zeros((nindex,ndays))
vtalktrue = np.zeros((nindex,ndays))

#NOTE Everything before this line can be computer ahead of time in our model
# Everything after this line is calculated for each room


################################################################################################
def room_pinf(ibinmax,kdecay,lowbreathe,maskfac_out,maskfac_in,
              roomvolx,achx,tendx,talkfractionx,vdotbreathex,
              vemtalk,vembreathe,kdep,fdep,
              nindex_room,nemmx_talk,nemmx_breathe,vinfx):

    # THIS ROUTINE CALCULATE THE PROBABILITY OF INFECTION IN A ROOM ASSUMING WELL-MIXED CONDITIONS AND FAR FIELD TRANSMISSION ONLY

    # NOTE
    # ibinmax,kdecay,lowbreathe,maskfac_out,maskfac_in are scalars that do not change from room to room
    # roomvolx,achx,tendx,talkfractionx,vdotbreathex are scalars that change from room to room 
    # vemtalk,vembreathe,kdep,fdep have dimension ibinmax and do no change from room to room
    # nemmx_talk,nemmx_breathe,vinfx have dimension nindex_room and change from room to room

    rnaconc = np.zeros((ibinmax,nindex_room))
    dose = 0.0

    for ibin in range(ibinmax):
        # Aerosol volume emission rate in each bin (mL/s)
        # nemmx_talk and nemmx_breathe are for lowbreathe  - scale by actual breathing
        vaerosol = talkfractionx * (nemmx_talk * vdotbreathex/lowbreathe) * vemtalk[ibin] + (1.0 - talkfractionx) * (nemmx_breathe * vdotbreathex/lowbreathe) * vembreathe[ibin] 

        # Loss coefficient (1/s)
        kloss = (kdecay + kdep[ibin] + achx)/3600.0 

        # CONCENTRATION IN WELL-MIXED CASE
        rnaconc[ibin,:] = ((1.0 - maskfac_out)*vinfx[:]*vaerosol)/(kloss*roomvolx) * (1.0 - (1 - np.exp(-kloss * tendx))/(kloss * tendx))

    # CALCULATE TOTAL RNA CONCENTRATION IN EACH BIN FROM ALL INDEX CASES IN ROOM
    totrnaconc = np.sum(rnaconc, axis=1)

    # CALCULATE DOSE FOR WELL MIXED CASE)
    for ibin in range(ibinmax):
        # DOSE FOR SUSCEPTIBLE IN WELL-MIXED CASE
        dose = dose + (1-maskfac_in)*fdep[ibin]*totrnaconc[ibin] * vdotbreathex * tendx

    # INFECTION PROBABILITY FOR SUSCEPTIBLE IN WELL MIXED CASE
    pinfx = 1.0 - np.exp(-dose/id50)

    return pinfx
################################################################################################


#NOTE Do not need anything after this line for our model
#Replace with network model

nrooms = nindex
for i_day in range(ndays):
    for i_room in range(nrooms):
        roomvolx = roomvol[i_room, i_day]
        achx = ach[i_room, i_day]
        tendx = tend[i_room, i_day]
        talkfractionx = talkfraction[i_room, i_day]
        vdotbreathex = vdotbreathe[i_room, i_day]
        nindex_room = 1
        nemmx_talk=np.zeros(nindex_room)
        nemmx_breathe=np.zeros(nindex_room)
        vinfx=np.zeros(nindex_room)
        for i_index in range(nindex_room):
            #assign a tag corresponding to the individual's index
            nemmx_talk[i_index] = nemm_talk[i_room] #Instead of i_room, it will be i_tag 
            nemmx_breathe[i_index] = nemm_breathe[i_room] #Same thing, use tag
            vinfx[i_index] = vinf_1[i_room,i_day] #replace iroom with tag
 
        pinf_1[i_room,i_day] = room_pinf(ibinmax,kdecay,lowbreathe,maskfac_out,maskfac_in,
              roomvolx,achx,tendx,talkfractionx,vdotbreathex,
              vemtalk,vembreathe,kdep,fdep,
              nindex_room,nemmx_talk,nemmx_breathe,vinfx)

        ninf_1[i_room,i_day] = np.rint(pinf_1[i_room,i_day] * nsus[i_room,i_day])



# CALCULATE R*
print('Number of new infections - WELL-MIXED; h=1.0 SALIVA: ', np.sum(ninf_1))
print('R*: ', np.sum(ninf_1)/nindex)
print('')

ninfindex_1 = np.sum(ninf_1, axis=1)

nzeros_1 = ninfindex_1[ninfindex_1 == 0]

ninfday_1 = np.sum(ninf_1, axis=0)

ezeros_1 = ninf_1.flatten()[ninf_1.flatten() == 0]

print('Max infections/index case WELL-MIXED; h=1.0 SALIVA: ',np.max(ninfindex_1))
print('')

print('Max infections/event WELL-MIXED; h=1.0 SALIVA: ',np.max(ninf_1))
print('')

print('Number of index cases with 0 infections WELL-MIXED; h=1.0 SALIVA: ', nzeros_1.size)
print('')

dzeros_1 = np.zeros(ndays)
print('Number of events with 0 infections WELL-MIXED; h=1.0 SALIVA: ', ezeros_1.size)
for iday in range(ndays):
    dzeros_1[iday] = ninf_1[:,iday][ninf_1[:,iday] == 0].size 
    print(iday,ninfday_1[iday],dzeros_1[iday])
print('')

