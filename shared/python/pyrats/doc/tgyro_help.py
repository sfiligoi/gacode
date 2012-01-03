"""This script gives help for the requested TGYROData data object."""

import sys
from pyrats.tgyro.data import TGYROData

if len(sys.argv) == 1:
    print "Please specify a directory containing TGYRO files."
else:
    sim1 = TGYROData(sys.argv[1])
    if len(sys.argv) == 2:
        print "Please specify a TGYROData method.  The options are:"
        print "data"
        print "docstrings (not a method, instead prints the docstrings for TGYROData)"
    else:
        if sys.argv[2] == 'data':
            for key in sorted(sim1.data.keys()):
                print key
        elif sys.argv[2] == 'docstrings':
            help(TGYROData)
        else:
            print "Incorrect TGYROData method.  The options are:"
            print "data"
            print "docstrings (not a method, instead prints the docstrings for TGYROData)"
