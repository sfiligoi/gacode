"""This script gives help for the requested profiles_genData data object."""

import sys
from pyrats.profiles_gen.data import profiles_genData

if len(sys.argv) == 1:
    print "Please specify a directory containing the files input.profiles and also possibly input.profiles.geo."
else:
    sim1 = profiles_genData(sys.argv[1])
    if len(sys.argv) == 2:
        print "Please specify a profiles_genData method.  The options are:"
        print "data"
        print "docstrings (not a method, instead prints the docstrings for profiles_genData)"
    else:
        if sys.argv[2] == 'data':
            print sorted(sim1.data.keys())
        elif sys.argv[2] == 'docstrings':
            help(profiles_genData)
        else:
            print "Incorrect profiles_genData method.  The options are:"
            print "data"
            print "docstrings (not a method, instead prints the docstrings for profiles_genData)"

