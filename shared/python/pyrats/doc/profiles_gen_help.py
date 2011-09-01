"""This script gives help for the requested profiles_genData data object."""

import sys
import pprint
from pyrats.profiles_gen.data import profiles_genData

if len(sys.argv) == 1:
    print "Please specify a directory containing input.profiles."
    sys.exit()

sim = profiles_genData(sys.argv[1])

if len(sys.argv) == 2:
    print "Please specify a profiles_genData method.  The options are:"
    print
    pprint.pprint(sorted(sim.__dict__.keys()))
else:
    if sys.argv[2] == 'data':
        for key in sorted(sim.data.keys()):
            print key.rjust(16),type(sim.data[key])
    elif sys.argv[2] == 'docstrings':
        help(profiles_genData)
    else:
        print "Incorrect profiles_genData method.  The options are:"
        print
        pprint.pprint(sorted(sim.__dict__.keys()))
        
