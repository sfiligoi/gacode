import sys
try:
    import matplotlib.pyplot as plt
except ImportError:
    print "This command requires matplotlib.  Please install matplotlib:"
    print "http://matplotlib.sourceforge.net/"
    sys.exit()
sys.path.append('/home/buuck/gacode/shared/python/profiles_gen')
from profiles_genData import profiles_genData
from math import *
from os.path import expanduser, expandvars

verbose = False
if int(sys.argv[5]) != 0:
    verbose = True
direc = sys.argv[1]
try:
    m1 = float(sys.argv[2])
except ValueError:
    print "ERROR: " + sys.argv[2] + " is not a valid number."
    sys.exit()
if m1 > 1:
    print u'ERROR: \u03c11 cannot be more than 1.'
    sys.exit()
if m1 < 0:
    print u'ERROR: \u03c11 cannot be less than 0.'
    sys.exit()
try:
    m2 = float(sys.argv[3])
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
    n = float(sys.argv[4])
except ValueError:
    print "ERROR: " + n + " is not a vaild number."
    sys.exit()
if n != floor(n):
    print "ERROR: Number of flux surfaces must be an integer."
    sys.exit()
min1 = min(m1, m2)
max1 = max(m1, m2)
if min1 == max1:
    print "ERROR: Radius interval must be greater than zero."
    sys.exit()
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
prof1.millerplot(min1, max1, n, verbose)
plt.show()
