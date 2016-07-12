import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

sim = cgyrodata('./')

ns = sim.n_species
t = sim.t

for moment in ['n','e']:
    if moment == 'n':
        mtag = '\Gamma'
        ttag = 'G'
        ftag = 'flux_n'
        if hasattr(sim,'kxky_flux_n'):
            y = np.sum(sim.kxky_flux_n,axis=(0,2))
        else:
            y = sim.flux_n
    elif moment == 'e':
        mtag = 'Q'
        ttag = 'Q'
        ftag = 'flux_e'
        if hasattr(sim,'kxky_flux_e'):
            y = np.sum(sim.kxky_flux_e,axis=(0,2))
        else:
            y = sim.flux_e

    data  = np.column_stack((sim.t,y[0,:]))
    head  = '(cs/a) t     '+ttag+'_1/'+ttag+'_GB'
    fname = 'out.cgyro.'+ftag
    for ispec in range(1,ns,1):
        head = head+'       '+ttag+'_'+str(ispec+1)+'/'+ttag+'_GB'
        data = np.column_stack((data,y[ispec,:]))
    np.savetxt(fname,data,fmt='%.8e',header=head)
    print 'INFO: (dump_flux_time) Created '+fname
