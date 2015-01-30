import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata


w = float(sys.argv[1])

sim = cgyrodata('./')

nt = sim.n_time

# Determine tmin
imin=0
for i in range(nt):
    if sim.t[i] < (1.0-w)*sim.t[-1]:
        imin = i+1

if imin == nt-1:
    print "Averaging Window too small." 
    sys.exit()

# Set print precision
np.set_printoptions(precision=3,suppress=False,threshold=100000)

print 'Average Window:',str(sim.t[imin])+' < (c_s/a) t < '+str(sim.t[-1])
print

b = np.zeros([sim.n_species])

for ispec in range(sim.n_species):
    s = 0.0
    for p in range(sim.n_n):
        s = s+average(sim.flux_n[ispec,p,:],sim.t,w)
    b[ispec] = s

print 'GAMMA [GB]  ',b
    
for ispec in range(sim.n_species):
    s = 0.0
    for p in range(sim.n_n):
        s = s+average(sim.flux_e[ispec,p,:],sim.t,w)
    b[ispec] = s

print 'Q     [GB]  ',b

