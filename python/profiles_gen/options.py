# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
import sys
from profiles_gen.data import profiles_genData

filevec = sys.argv[1].split(',')

prof = profiles_genData(filevec[0])
keys = sorted(prof.data.keys())

print() 
for i in range(len(keys)):
    print(keys[i].split()[0])

keys = sorted(prof.geo.keys())

print() 
for i in range(len(keys)):
    print(keys[i].split()[0])
