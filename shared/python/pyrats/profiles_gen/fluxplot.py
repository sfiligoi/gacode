import sys
import matplotlib.pyplot as plt
from math import *
from os.path import expanduser, expandvars
from pyrats.profiles_gen.data import profiles_genData

verbose = False

simdir = sys.argv[1]

try:
    m1 = float(sys.argv[2])
except ValueError:
    print "ERROR: " + sys.argv[2] + " is not a valid number."
    sys.exit()

if m1 > 1.0 or m1 < 0.0:
    print u'ERROR: \u03c1 out of bounds.'
    sys.exit()

try:
    m2 = float(sys.argv[3])
except ValueError:
    print "ERROR: " + m2 + " is not a valid number."
    sys.exit()

if m2 > 1.0 or m2 < 0.0:
    print u'ERROR: \u03c1 out of bounds.'
    sys.exit()

try:
    n = int(sys.argv[4])
except ValueError:
    print "ERROR: " + n + " is not a valid number."
    sys.exit()

minr = min(m1, m2)
maxr = max(m1, m2)

if minr == maxr:
    print "ERROR: Radius interval must be greater than zero."
    sys.exit()

try:
    prof = profiles_genData(simdir)
except IOError:
    print "ERROR: No profile data in "+simdir 
    sys.exit()

maxn = prof.match(maxr,prof.get('rho (-)'))-prof.match(minr,prof.get('rho (-)'))

if maxn == 0:
    print "ERROR: No flux surfaces in given interval."
    sys.exit()

if maxn < n:
    print "Warning: Maximum number of flux surfaces for this interval"
    print "is " + str(maxn) + ".  Changing n to " + str(maxn) + "."
    n = maxn

if int(sys.argv[5]) == 1:
    prof.millerplot(minr, maxr, n, verbose)
if int(sys.argv[5]) == 2:
    prof.fourierplot(minr, maxr, n, verbose)
if int(sys.argv[5]) == 3:
    prof.compplot(minr, maxr, n, verbose)

plt.show()
