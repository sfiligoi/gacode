# Dump text report of flux

import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata


w = float(sys.argv[1])
ext = sys.argv[2]

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

print 'INFO: Average Window:',str(sim.t[imin])+' < (c_s/a) t < '+str(sim.t[-1])

if ext == 'dump':
    # Datafile output

    b = np.zeros([2,sim.n_n])

    fn=open('out.cgyroplot.fluxn','w')
    fe=open('out.cgyroplot.fluxe','w')
    for ispec in range(sim.n_species):
        fn.write(str(ispec)+' ')
        fe.write(str(ispec)+' ')
        for i_n in range(sim.n_n):
            yn = np.sum(sim.flux_n,axis=0)
            ye = np.sum(sim.flux_e,axis=0)
            b[0,i_n] = average(yn[ispec,i_n,:],sim.t,w)
            b[1,i_n] = average(ye[ispec,i_n,:],sim.t,w)
            fn.write("{:.3e}".format(b[0,i_n])+' ')
            fe.write("{:.3e}".format(b[1,i_n])+' ')
        fn.write('\n')
        fe.write('\n')            
        print 'Wrote output to out.cgyroplot.flux*'
else:
    # Screen output

    # Set print precision
    np.set_printoptions(precision=3,suppress=False,threshold=100000)

    b = np.zeros([sim.n_species])

    for ispec in range(sim.n_species):
        y = np.sum(sim.flux_n,axis=(0,2))
        b[ispec] = average(y[ispec,:],sim.t,w)

    print 'GAMMA [GB]  ',b
    
    for ispec in range(sim.n_species):
        y = np.sum(sim.flux_e,axis=(0,2))
        b[ispec] = average(y[ispec,:],sim.t,w)

    print 'Q     [GB]  ',b


    
