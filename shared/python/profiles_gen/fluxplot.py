"""fluxplot.py contains plotting routines for data parsed with
profiles_genData.  As the name suggests, it plots flux surfaces.

This module mostly just checks for errors in the user input and then passes
control over to profiles_genData.

The module requires at least 4 command line arguments to run: the module
itself, a directory location, a plot type, and a radius to plot.  Optionally,
the user can specify another radius, and a number of radii to plot between
the two given radii.  See ~/gacode/shared/bin/profiles_gen for more
information on the command line interface."""

import sys
import matplotlib.pyplot as plt
from profiles_genData import profiles_genData
from math import *
from os.path import expanduser, expandvars

#Default values
direc = '.'
typ = '-c'
min1 = 1
max1 = 0
n = 1
l = 4
prof1 = 0

#Error catching
if len(sys.argv) > 1:
    direc = sys.argv[1]
if len(sys.argv) > 2:
    typ = sys.argv[2]
if len(sys.argv) > 3:
    try:
        m1 = float(sys.argv[3])
    except ValueError:
        print "ERROR: " + sys.argv[3] + " is not a valid number."
        sys.exit()
    if m1 > 1:
        print u'ERROR: \u03c11 cannot be more than 1.'
        sys.exit()
    if m1 < 0:
        print u'ERROR: \u03c11 cannot be less than 0.'
        sys.exit()

#If there are 6 arguments, then the user is looking for a range of radii to be
#plotted.
if len(sys.argv) == 6:
    #More error catching
    try:
        m2 = float(sys.argv[4])
    except ValueError:
        print "ERROR: " + m2 + " is not a vaild number."
        sys.exit()
    if m2 > 1:
        print u'ERROR: \u03c12 cannot be more than 1.'
        sys.exit()
    if m2 < 0:
        print u'ERROR: \u03c12 cannot be less than 0.'
        sys.exit()
    try:
        n = float(sys.argv[5])
    except ValueError:
        print "ERROR: " + n + " is not a vaild number."
        sys.exit()
    if n != floor(n):
        print "ERROR: Number of flux surfaces must be an integer."
        sys.exit()
    if n < 2:
        print "ERROR: Number of flux surfaces must at least 2."
        sys.exit()
    min1 = min(m1, m2)
    max1 = max(m1, m2)
    if min1 == max1:
        print "ERROR: Radius interval must be greater than zero."
        sys.exit()
    #Checks to see if requested number of surfaces is greater than total
    #number available.
    try:
        prof1 = profiles_genData(direc)
    except IOError:
        print "ERROR: " + expanduser(expandvars(direc))
        print "does not contain file input.profiles and/or file input.profiles.geo.  Type profiles_gen for help."
        sys.exit()
    maxn = prof1.match(max1, prof1.get('rho (-)')) - prof1.match(min1, prof1.get('rho (-)'))
    if maxn == 0:
        print "ERROR: No flux surfaces in given interval.  Please specify a larger interval."
        sys.exit()
    if maxn < n:
        print "Warning: Maximum number of flux surfaces for this interval"
        print "is " + str(maxn) + ".  Changing n to " + str(maxn) + "."
        n = maxn
    l = 6

#If there are only 4 arguments, then the user is only looking for one radius.
if len(sys.argv) == 4:
    min1 = m1

#Any other number of arguments is an error.
if len(sys.argv) == 5 or len(sys.argv) > 6:
    print "ERROR: Wrong number of arguments.  Type profiles_gen for help."
    sys.exit()

if not prof1:
    try:
        prof1 = profiles_genData(direc)
    except IOError:
        print "ERROR: " + expanduser(expandvars(direc))
        print "does not contain file input.profiles and/or file input.profiles.geo.  Type profiles_gen for help."
        sys.exit()
prof1.fplot(typ, l, min1, max1, n)
plt.show()
