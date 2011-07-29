"""plot.py contains plotting routines for data parsed with
profiles_genData.

This module mostly checks for errors, and contains code to print all possible
variables to the terminal.

The default values for this module are: directory='.', parameters=ne, ni_1, Te, 
Ti_1, dimensions=2x2.  Optionally, the user can specify an arbitrary number of
variables to plot, as well as the dimensions of the plotting windows.  If -c
is entered as the second argument (after the name of the module itself), then
the user can specify two directories for comparison.  See
~/gacode/shared/bin/profiles_gen for more information on the command line
interface."""

import sys
from pyrats.profiles_gen.data import profiles_genData
from errorcheck import *
try:
    import matplotlib.pyplot as plt
except ImportError:
    print "This command requires matplotlib.  Please install matplotlib:"
    print "http://matplotlib.sourceforge.net/"
    sys.exit()
import math

args = ['ne', 'ni_1', 'Te', 'Ti_1']
prof1 = opendir(sys.argv[1])

n1 = setdim(sys.argv[2])
n2 = setdim(sys.argv[3])

if len(sys.argv) < 5:
    for arg in args:
        checkarg(prof1, arg)
        prof1.plot(arg, n1, n2)
    plt.show()
else:
    args = sys.argv[4:]
    for arg in args:
        checkarg(prof1, arg)
        prof1.plot(arg, n1, n2)
    plt.show()
