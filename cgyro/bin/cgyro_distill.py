#!/usr/bin/env python

import sys
import os
import numpy as np
from gacodefuncs import *
from cgyro.data import cgyrodata

if len(sys.argv) < 2:
    print "Usage: cgyro_distill.py <window>"
    sys.exit()

if not os.path.isfile('bin.cgyro.kxky_phi'):
    print "Error: no data found"
    sys.exit()
    
w = float(sys.argv[1])

# read basic simulation data
sim = cgyrodata('./')

# read larger files
sim.getbigfield()

nt  = sim.n_time
nr  = sim.n_radial
nn  = sim.n_n
ns  = sim.n_species
nth = sim.theta_plot

t = sim.t

# Get index for average window
imin,imax=iwindow(t,w,0.0)

header = 'Averaging over '+str(t[imin])+' < (c_s/a) t < '+str(t[imax])
header = header+' nkx = '+str(nr)+' nky = '+str(nn)+'\n'

print(header)

# kxky_phi -> [2,self.n_radial,self.theta_plot,self.n_n,nt]

# 1. Average phi^2
phi2 = sim.kxky_phi[0,:,0,:,:]**2+sim.kxky_phi[1,:,0,:,:]**2

ave_phi2 = np.zeros([nr,nn])
for j in range(nn):
    ave_phi2[:,j] = average_n(phi2[:,j,:],t,w,0.0,nr)

np.savetxt('out.cgyro.phi2kxky',ave_phi2,fmt='%1.6e',header=header)

print('Wrote out.cgyro.phi2kxky')

# 2. Average Re(phi_0)
phi0 = sim.kxky_phi[:,:,0,0,:]
ave_phi0 = np.zeros([nr,2])
ave_phi0[:,0] = average_n(phi0[0,:,:],t,w,0.0,nr)
ave_phi0[:,1] = average_n(phi0[1,:,:],t,w,0.0,nr)

np.savetxt('out.cgyro.phi0kx',ave_phi0,fmt='%1.6e',header=header)

print('Wrote out.cgyro.phi0kx')
