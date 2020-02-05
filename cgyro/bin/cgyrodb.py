#!/usr/bin/env python

import os
#import yaml
import numpy as np
import pandas as pd
from gacodeinput import *
from gacodefuncs import *
from cgyro.data import cgyrodata

mode = 'flux'

meta = {}

def gendict(sim,mode):

   mydict = {}

   mydict['n_radial'] = sim.n_radial
   mydict['n_theta'] = sim.n_theta

   if mode == 'flux':
      sim.getflux()
      y = np.sum(sim.ky_flux,axis=(2,3))   
      for ispec in range(sim.n_species):
         g = average(y[ispec,0,:],sim.t,0.5,0.0)
         q = average(y[ispec,1,:],sim.t,0.5,0.0)
         v = average(y[ispec,2,:],sim.t,0.5,0.0)
         mydict['g'+str(ispec)] = g
         mydict['q'+str(ispec)] = q
         mydict['v'+str(ispec)] = v

   return mydict

for mdir in next(os.walk('./'))[1]:
   sim = cgyrodata(mdir+'/',silent=True)
   meta[mdir]=gendict(sim,mode)

print(meta)

df = pd.DataFrame(meta) 
print(df.T[['n_radial','q0']])

#with open('data.yml','w') as outfile:
#   yaml.dump(meta,outfile,default_flow_style=False)
