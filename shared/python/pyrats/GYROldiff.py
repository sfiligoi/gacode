from pyrats.data import GYROData
import sys

verbose = bool(int(sys.argv[2]))
sim1 = GYROData(sys.argv[1])

sim1.make_diff()
if verbose:
    print "Particle and energy diffusivity:"
    print sim1.diff
else:
    print "Energy diffusivity for kinetic species zero, and field zero:"
    print sim1.diff['chi_sigma/chi_{GB}'][0, 0, :]
