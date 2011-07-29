from pyrats.gyro.data import GYROData
import sys
import numpy as np

verbose = bool(int(sys.argv[2]))
sim1 = GYROData(sys.argv[1])
n = int(sys.argv[3])-1

sim1.make_gbflux()
temp = np.sum(sim1.gbflux, axis=1)

if n > int(sim1.profile['n_spec']):
    print "Warning: Max number of species is " + str(int(sim1.profile['n_spec'])) + ".",
    print " Changing number of species to " + str(int(sim1.profile['n_spec']))
    n = int(sim1.profile['n_spec']) - 1

if verbose:
    print
    print "Gyrobohm-normalized particle and energy flux:"
    print
    print "      TIME  ",
    for a in range(n):
        print "   |   PARTICLE FLUX SPE " + str(a) + "   |    ENERGY FLUX SPE " + str(a),
    print "   |   PARTICLE FLUX SPE " + str(n) + "   |    ENERGY FLUX SPE " + str(n)    
    for i in range(len(temp[0, 0, :])):
        print repr(i).rjust(15), '|',
        for a in range(n):
            print repr(temp[a, 0, i]).ljust(23), '|', repr(temp[a, 1, i]).ljust(23), '|',
        print repr(temp[n, 0, i]).ljust(23), '|', repr(temp[n, 1, i]).ljust(23)
else:
    print
    print "Gyrobohm-normalized energy flux for species zero:"
    print
    print "      TIME      |       ENERGY FLUX"
    for i in range(len(temp[0, 0, :])):
        print repr(i).rjust(15), '|', repr(temp[0, 1, i]).ljust(23)
