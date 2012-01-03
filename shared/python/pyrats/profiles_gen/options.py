import sys
import string
from pyrats.profiles_gen.data import profiles_genData

prof = profiles_genData(sys.argv[1])
keys = sorted(prof.data.keys())

print 
for i in range(len(keys)):
    print keys[i].split()[0]
