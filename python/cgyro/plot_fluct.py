import sys
import numpy as np
from gacodeplotdefs import *
from gacodefuncs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
moment = sys.argv[2]
ymin = sys.argv[3]
ymax = sys.argv[4]

sim = cgyrodata('./')
nt = sim.n_time

#(2,n_radial,n_species,n_n,nt)
n_chunk = 2*sim.n_radial*sim.n_species*sim.n_n*nt
for line in open('out.cgyro.kxky_n'):
    print float(line)

print "INFO: (data.py) Read data"
