"""profiles_generrorcheck.py contains routines that are useful to other python scripts to catch
errors in input.""" 

import sys
from pyrats.profiles_gen.data import profiles_genData
import matplotlib.pyplot as plt
import math

#First are several functions that perform error checking that comes up frequently.
def opendir(directory):
    try:
        prof = profiles_genData(directory)
    except IOError:
        print directory + " does not contain file input.profiles.  Type profiles_gen for help."
        sys.exit()
    except IndexError:
        print "ERROR: Too few arguments.  Type profiles_gen for help."
        sys.exit()
    return prof

def setdim(n):
    try:
        rval = int(n)
    except ValueError:
        print "ERROR: " + n + " is not a valid dimension.  Type profiles_gen for help."
        sys.exit()
    except IndexError:
        print "ERROR: Please give second dimension.  Type profiles_gen for help."
        sys.exit()
    return rval

def checkarg(prof, arg):
    flag = 0
    for k in prof.data.iterkeys():
        if arg == k.split()[0]:
            flag = 1
    if not flag:
        print "ERROR: ", arg, " is not a valid parameter.  Type -options after -plot for help."
        sys.exit()

#getoptions is not error checking, but is still a routine that gets used
#frequently.
def getoptions(prof):
    keys = []
    for k, v in prof.data.iteritems():
        s = 0
        for item in v:
            s = s + float(item)
        if s != 0:
        #Unless every entry in v is zero, it gets added to the list.
            keys.append(k.split()[0])
    keys.sort()
    for k in keys:
        print k
