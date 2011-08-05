"""This script gives help for the requested NEOData data object."""

import sys
from pyrats.neo.data import NEOData

if len(sys.argv) == 1:
    print "Please specify a directory containing NEO files."
else:
    sim1 = NEOData(sys.argv[1])
    if len(sys.argv) == 2:
        print "Please specify a NEOData method.  The options are:"
        print "transport"
        print "HH_theory"
        print "CH_theory"
        print "TG_theory"
        print "S_theory"
        print "HR_theory"
        print "HS_theory"
        print "docstrings (not a method, instead prints the docstrings for NEOData)"
    else:
        if sys.argv[2] == 'transport':
            for key in sorted(sim1.transport.keys()):
                print key
        elif sys.argv[2] == 'HH_theory':
            for key in sorted(sim1.HH_theory.keys()):
                print key
        elif sys.argv[2] == 'CH_theory':
            for key in sorted(sim1.CH_theory.keys()):
                print key
        elif sys.argv[2] == 'TG_theory':
            for key in sorted(sim1.TG_theory.keys()):
                print key
        elif sys.argv[2] == 'S_theory':
            for key in sorted(sim1.S_theory.keys()):
                print key
        elif sys.argv[2] == 'HR_theory':
            for key in sorted(sim1.HR_theory.keys()):
                print key
        elif sys.argv[2] == 'HS_theory':
            for key in sorted(sim1.HS_theory.keys()):
                print key
        elif sys.argv[2] == 'docstrings':
            help(NEOData)
        else:
            print "Incorrect NEOData method.  The options are:"
            print "transport"
            print "HH_theory"
            print "CH_theory"
            print "TG_theory"
            print "S_theory"
            print "HR_theory"
            print "HS_theory"
            print "docstrings (not a method, instead prints the docstrings for NEOData)"
