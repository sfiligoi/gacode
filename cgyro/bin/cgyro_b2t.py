# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
import sys
import numpy as np

if len(sys.argv) == 1:
    print('python cgyro_b2t.py <binary-source> <text-dest>')
    sys.exit()
    
bfile = sys.argv[1]
afile = sys.argv[2]

arr = np.fromfile(bfile,dtype=np.float64)

print("Read binary data in "+bfile)

np.savetxt(afile,arr,fmt='%.5E')

print("Wrote ASCII data to "+afile)
