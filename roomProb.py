import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as colors
import random

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

rng = default_rng(190319381062290118519841602886838153481)
float_epsilon = np.finfo(float).eps
current_day = 0 #Will make this a parameter when we turn the entire thing into a function

#Map each person to a viral load profile
#Do this in external network model
vloads = {}
for person in range(len(model)):
    idx = random.randint(0,10001)
    vloads[person] = idx



#s = 0 
#e = 1
#i = 2
#r = 3

#Model is an array with people, status, and days since infected
n = 100000 #100,000 people
model = np.zeros((n, 2))

#Going to make a test case here
for i in range(100000):
    num = random.randint(1,100000)
    if num <= 10000:
        model[i][0] = 2
        model[i][1] = 0
    else:
        model[i][0] = 0
        model[i][1] = float("inf")



room = random.sample(range(1,100000), 100) #person indices
#from this we can get number of index cases in the room, breathe and talk rates, and viral load (since each person will be mapped to a viral load)
#all other variables are calculated externally 
nindex_room = len([num for num in room if model[num][0] == 2])
areapp = random.uniform(1,10)
roomarea = (1 + (len(room)-nindex_room))*areapp
roomheight = 3.0
roomvol = roomheight*roomarea

meanach = 2.00
stdach = 2.45

mu_ach = np.log(meanach*meanach/np.sqrt(meanach*meanach + stdach*stdach))
sigma_ach = np.sqrt(np.log(1.0 + stdach*stdach/(meanach*meanach)))
ach = rng.lognormal(mu_ach,sigma_ach)

mu_tend = np.log(1.0)
sigma_tend=0.5
lntend_min = np.log(0.1)
lntend_max = np.log(3.0+float_epsilon)
aclip = (lntend_min - mu_tend)/sigma_tend
bclip = (lntend_max - mu_tend)/sigma_tend

lntend_clip = truncnorm.rvs(aclip, bclip, loc=mu_tend, scale=sigma_tend,random_state=rng)
tend = 3600.0*np.exp(lntend_clip)

talkfraction = rng.uniform(low=0.05,high=0.25+float_epsilon)


lowbreathe = 0.51/3600.0
highbreathe = 1.24/3600.0
# lowbreathe = 1.0e-4
# highbreathe = 2.0e-4

# Implicit assumption here is that breathing rate is same for all individuals in room
br_activity = [0.51/3600.0, 0.57/3600.0, 1.24/3600]

activityindex = rng.choice(3)
vdotbreathe=br_activity[activityindex]

mean = 3.4*13
std = 2.7*13
mu = np.log(mean*mean/np.sqrt(mean*mean + std*std))
sigma = np.sqrt(np.log(1.0 + std*std/(mean*mean)))
lmax = np.log(20.0*13)
bclip = (lmax - mu)/sigma
lclip = truncnorm.rvs(-np.inf, bclip, loc=mu, scale=sigma, size=len(model),random_state = rng)
nemm_talk = np.exp(lclip)
nemm_talk_room = [nemm_talk[i] for i in room if model[i][0] == 2]

mub = np.log(8) - 0.48*np.log(10)
sigmab = 0.85*np.log(10)/(np.sqrt(2))
nemm_breathe = rng.lognormal(mub, sigmab, size=len(model))
nemm_breathe_room = [nemm_breathe[i] for i in room if model[i][0] == 2]

kevi = np.loadtxt('data/Ke_Infectiousness_10000.csv', delimiter = ',', skiprows = 1)

#Build viral load list


vinf = []
for person in room:
    if model[person][0] == 2:
        load = kevi[int(current_day - model[person][1])][vloads[person]]
        vinf.append(load)


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
    #nindex_room is the number of index cases in the room
    #can fill this value at room assignment

    dose = 0.0

    for ibin in range(ibinmax):
        # Aerosol volume emission rate in each bin (mL/s)
        # nemmx_talk and nemmx_breathe are for lowbreathe  - scale by actual breathing
        vaerosol = talkfractionx * (nemmx_talk * vdotbreathex/lowbreathe) * vemtalk[ibin] + (1.0 - talkfractionx) * (nemmx_breathe * vdotbreathex/lowbreathe) * vembreathe[ibin] 
        #talk fraction comes from uniform -> can be assigned ahead of time
        #
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