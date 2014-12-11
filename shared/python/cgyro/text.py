import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

data = cgyrodata(sys.argv[1])
w    = float(sys.argv[2])

s = 0.0
for p in range(data.n_n):
    f = average(data.flux_e[p,:],data.t,w)
    s = s+f
    print p,f
    
print 'TOTAL: ',s
