#!/usr/bin/python
#For calculating timings required in multiple resonance SHARPER experiments with two signals
#Suitable for two off-resonance peaks with different |Frequency|
#Calculates a timing for a spectral width and refocussing period within user given bound
#Bounds set to provide nice spectra with adequate decoupling and T2* - > T2
#Wider bounds will improve average spectral quality
#N.B. FutureWarnings are harmless, but may render this way of performing the calcuations obsolete in the future
##N.B. This script calculates the value of "Delta" for multiple resonance SHARPER which is the half chunk duration, not the full chunk duration
import math
import numpy as np
import warnings
from scipy import special
warnings.filterwarnings("ignore", category=FutureWarning) #This ignores a depreciation warning that is currently irrelevant

#####Start of parameters the user should change#####
#Frequency of NMR signal (Hz)
f1 = 1907.9690
#Spectrometer operating Frequency (for the nucleus of interest!)
Spectrometer_Frequency = 500.0 #This is in MHz, for the nucleus of interest
#Allowable Spectral Width
ppm_upper = 40.0
ppm_lower = 10.0
#Allowable chunk timings in milliseconds (enter as a decimal, needs to be a float)
longest_timing = 8.0
shortest_timing = 4.0
#####End of parameters the user should change#####

#Calculate the period of the NMR signal
p1 = 1.0/f1
#Convert to Hz limits
hz_upper=(ppm_upper*Spectrometer_Frequency)
hz_lower=(ppm_lower*Spectrometer_Frequency)
#Calculate the lowest allowed dwell time
dwell_lower=(1e6*(0.5/hz_upper))
#Calculate the maximum allowed dwell time
dwell_upper=(1e6*(0.5/hz_lower))
#Fetch Array Files
t_ar = np.genfromtxt("TimingArray.csv",delimiter=",")
dwell_ar = np.genfromtxt("Dwellarray.csv",delimiter=",")
cp_ar = np.genfromtxt("ChunkPointarray.csv",delimiter=",")

#Produces a boolean mask for suitable dwell times
b = ((dwell_upper>dwell_ar)&(dwell_ar>dwell_lower))
#Produce arrays modified by first boolean mask with only suitable dwell times
dwell_arb = dwell_ar[b]
cp_arb = cp_ar[b]
t_arb = t_ar[b]
#Produce a second boolean mask based on timings
c = (((shortest_timing/1000)<t_arb)&(t_arb<(longest_timing/1000)))
#Produce arrays modified by second boolean mask with only suitable dwell times
dwell_arbc = dwell_arb[c]
cp_arbc = cp_arb[c]
t_arbc = t_arb[c]

#runs a modulo operation on the timing array with boolean mask applied
mod_arf1 = np.mod(t_arbc, p1) 
#turn into a spectral quality array
mod_sqp1 = np.divide(mod_arf1, p1/100)
#Modify arrays to allow incomplete rotations
sub50ind = [(50.00<mod_sqp1)*(mod_sqp1<75.00)]
mod_sqp1[sub50ind] -= 50.0
subfrom50ind = [(25.00<mod_sqp1)*(mod_sqp1<50.00)]
mod_sqp1[subfrom50ind] -= 50.0 
subfrom100ind = [(75.00<mod_sqp1)*(mod_sqp1<100.00)]
mod_sqp1[subfrom100ind] -= 100.0 
mod_sqp1 = np.absolute(mod_sqp1)
#Troubleshooting line to print arrays
#np.savetxt('fname.txt',mod_sqp1)
#Find minimum element in modulo array
minElement = np.amin(mod_sqp1)
index = np.where(mod_sqp1 == np.amin(mod_sqp1))
#N.B. In the case multiple experimental viable timings are equally close to theoretically valid timings it will report them all
#Fetch the corresponding Chunkpoints and dwell time
print "OptimalDwellTime =", dwell_arbc[index]
print "ChunkpointsMultiple =", cp_arbc[index]
#print index
sq1 = mod_sqp1[index]
#Spectral Quality Prediction - Percentage of final rotation completed
#If this number is not close to 0 or 100 the spectra can be expected to be of poor quality (S/A ratio)
Spectral_Quality_f1 = (sq1)
print 'Percentage of final rotation completed =', Spectral_Quality_f1,'%'
