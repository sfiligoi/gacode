from pyrats.data import GYROData
import sys

verbose = bool(int(sys.argv[2]))
sim1 = GYROData(sys.argv[1])

sim1.make_gbflux()
if verbose:
    print "Gyro-Bohm normalized particle, energy, and momentum flux, and exchange power density:"
    print sim1.gbflux
else:
    print "Gyro-Bohm normalized energy flux for species zero, field zero:"
    print sim1.gbflux[0, 0, 1, :]
