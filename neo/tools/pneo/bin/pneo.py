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
#  ingeodata(1)  = I/psiprime
#  ingeodata(2)  = ftrap
#  ingeodata(3)  = <B^2>
#  ingeodata(4)  = <1/B^2>
#  ingeodata(5)  = <(b dot grad B)^2>
#  ingeodata(6)  = <(b dot grad B)^2/B^2>
#  ingeodata(7)  = <q R0 / Jpsi B>
#  ingeodata(8)  = <Btor^2>
#  ingeodata(9)  = <Bpol^2>
#  ingeodata(10) = <|grad r|^2>
#  ingeodata(11) = <|grad r|^2/B^2>
#  ingeodata(12) = <-grad_r * gsin/B>

ingeodata = np.fromfile('ingeodata.dat',dtype='float',sep=" ")
n = len(ingeodata)/12
ingeodata = np.reshape(ingeodata,(12,n),'F')

# Outputs
#  outdata(1) = Cne
#  outdata(2) = CTe
#  outdata(3) = Cni1
#  outdata(4) = CTi1
#  outdata(5) = Cni2
#  outdata(6) = CTi2

outdata = np.fromfile('outdata.dat',dtype='float',sep=" ")
n = len(indata)/6
outdata = np.reshape(outdata,(6,n),'F')


