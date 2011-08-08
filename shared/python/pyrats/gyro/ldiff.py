"""This file is executed by the bash script gyro_plot when a listing of the
diffusivities is requested."""

from pyrats.gyro.data import GYROData
import sys
import numpy as np

verbose = bool(int(sys.argv[2]))
sim1 = GYROData(sys.argv[1])
n = int(sys.argv[3])-1

sim1.make_diff()

if n > int(sim1.profile['n_spec']):
    print "Warning: Max number of species is " + str(int(sim1.profile['n_spec'])) + ".",
    print " Changing number of species to " + str(int(sim1.profile['n_spec']))
    n = int(sim1.profile['n_spec']) - 1

if verbose:
    print
    print "Gyrobohm-normalized particle and energy diffusivities averaged over radius and summed over mode number:"
    print
    print "      TIME      ",
    for a in range(n):
        print "|PARTICLE DIFFUSIVITY SPE " + str(a) + "|ENERGY DIFFUSIVITY SPE " + str(a),
    print "|PARTICLE DIFFUSIVITY SPE " + str(n) + "|ENERGY DIFFUSIVITY SPE " + str(n)
    temp = []
    temp = np.sum(sim1.diff, axis=1)
    for i in range(len(temp[0][0])):
        print repr(i).rjust(16), '|',
        for a in range(n):
            print repr(temp[0][a][i]).ljust(24), '|', repr(temp[1][a][i]).ljust(23), '|',
        print repr(temp[0][n][i]).ljust(24), '|', repr(temp[1][n][i]).ljust(23)
else:
    print
    print "Gyrobohm-normalized particle and energy diffusivities averaged over radius and summed over mode number for species 0:"
    print
    print "      TIME       |   ENERGY DIFFUSIVITY"
    temp = []
    temp = np.sum(sim1.diff, axis=1)
    for i in range(len(temp[0][0])):
        print repr(i).rjust(16), '|', repr(temp[0][1][i]).ljust(23)
