import os
import numpy as np
import sys

# Inputs
#  indata(1) = neo_rmin_over_a_in
#  indata(2) = neo_q_in
#  indata(3) = neo_nu_1_in
#  indata(4) = neo_dens_in(2)
#  indata(5) = neo_temp_in(2)
#  indata(6) = neo_delta_in
#  indata(7) = neo_s_delta_in
#  indata(8) = neo_s_kappa_in

indata = np.fromfile('indata.dat',dtype='float',sep=" ")
n = len(indata)/8
indata = np.reshape(indata,(8,n),'F')

# Geometry inputs
#  ingeodata(1) = 

ingeodata = np.fromfile('ingeodata.dat',dtype='float',sep=" ")
n = len(ingeodata)/12
ingeodata = np.reshape(ingeodata,(12,n),'F')

# Outputs
#  outdata(1) = 

outdata = np.fromfile('outdata.dat',dtype='float',sep=" ")
n = len(indata)/6
outdata = np.reshape(outdata,(6,n),'F')


