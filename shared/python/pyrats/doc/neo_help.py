"""This script gives help for the requested NEOData data object."""

import sys
import pprint
import numpy as np
from pyrats.neo.data import NEOData

if len(sys.argv) == 1:
    print "Please specify a directory containing NEO files."
    sys.exit()

sim = NEOData(sys.argv[1])
if len(sys.argv) == 2:
    print "Please specify a NEOData method.  The options are:"
    print
    pprint.pprint(sorted(sim.__dict__.keys()))
else:
    if sys.argv[2] == 'dirname':
        print 'dirname',type(sim.dirname)
    elif sys.argv[2] == 'equil':
        for key in sorted(sim.equil.keys()):
            print key.rjust(16),type(sim.equil[key])
    elif sys.argv[2] == 'grid':
        for key in sorted(sim.grid.keys()):
            print key.rjust(16),type(sim.grid[key])
    elif sys.argv[2] == 'rotation':
        for key in sorted(sim.rotation.keys()):
            print key.rjust(16),type(sim.rotation[key])
    elif sys.argv[2] == 'theory':
        for key in sorted(sim.theory.keys()):
            print key.rjust(16),type(sim.theory[key])
    elif sys.argv[2] == 'transport':
        for key in sorted(sim.transport.keys()):
            print key.rjust(16),type(sim.transport[key])
    elif sys.argv[2] == 'transport_exp':
        for key in sorted(sim.transport_exp.keys()):
            print key.rjust(16),type(sim.transport_exp[key])
    elif sys.argv[2] == 'vel':
        print 'vel:',sim.veltag,np.shape(sim.vel)
    elif sys.argv[2] == 'docstrings':
        help(NEOData)
    else:
        print "Incorrect NEOData method.  The options are:"
        print
        pprint.pprint(sorted(sim.__dict__.keys()))
