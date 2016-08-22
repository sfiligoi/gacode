# Text-based checker for TGYRO data class

import sys
import numpy as np
from gacodeplotdefs import *
from tgyro.ndata import tgyrodata

sim = tgyrodata(sys.argv[1])

print 'Number of ions  : ',sim.n_ion
print 'Number of radii : ',sim.n_r
print 'Evolution eqns  : ',sim.n_evolve
print 'Completed iter  : ',sim.n_iterations

p=0
for key in sorted(sim.data.keys()):
    p=p+1
    print p,' : ',key
