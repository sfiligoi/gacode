import sys
import string
from pyrats.profiles_gen.data import profiles_genData

filevec = string.splitfields(sys.argv[1],',')

prof = profiles_genData(filevec[0])
keys = sorted(prof.data.keys())

print 
for i in range(len(keys)):
    print keys[i].split()[0]

keys = sorted(prof.geo.keys())

print 
for i in range(len(keys)):
    print keys[i].split()[0]
