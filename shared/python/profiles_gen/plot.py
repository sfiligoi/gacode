"""plot.py contains plotting routines for data parsed with
profiles_genData.

This module mostly checks for errors, and contains code to print all possible
variables to the terminal.

The module requires at least 3 command line arguments to run: the module itself,
a directory location, and a variable to plot.  Optionally, the user can specify
an arbitrary number of variables to plot, as well as the dimensions of the
plotting windows.  See ~/gacode/shared/bin/profiles_gen for more information on
the command line interface."""

import sys
from profiles_genData import profiles_genData
import matplotlib.pyplot as plt

args = []
n1 = 1
n2 = 1

#Error catching
try:
    prof1 = profiles_genData(sys.argv[1])
except IOError:
    print sys.argv[1]
    print "does not contain file input.profiles.  Type profiles_gen for help."
    sys.exit()
except IndexError:
    print "ERROR: Too few arguments.  Type profiles_gen for help."
    sys.exit()

if len(sys.argv) > 2:
    #This code segment prints all of the variables with nonzero entries when
    #'-options' is called.
    if sys.argv[2] == '-options':
        keys = []
        for k, v in prof1.data.iteritems():
            s = 0
            for item in v:
                s = s + float(item)
            if s != 0:
                #Unless every entry in v is zero, it gets added to the list.
                keys.append(k.split()[0])
        keys.sort()
        for k in keys:
            print k
    #More error catching
    elif sys.argv[2].isdigit():
        args = sys.argv[4:]
        try:
            n1 = int(sys.argv[2])
        except ValueError:
            print "ERROR: " + sys.argv[2] + " is not a valid dimension.  Type profiles_gen for help."
            sys.exit()
        try:
            n2 = int(sys.argv[3])
        except ValueError:
            print "ERROR: " + sys.argv[3] + " is not a valid dimension.  Type profiles_gen for help."
            sys.exit()
        for arg in args:
            flag = 0
            for k in prof1.data.iterkeys():
                if arg == k.split()[0]:
                    flag = 1
            if not flag:
                print "ERROR: ", arg, " is not a valid parameter.  Type -options after -plot for help."
                sys.exit()
            prof1.plot(arg, n1, n2)
        plt.show()
    elif sys.argv[2].isdigit() == False:
        args = sys.argv[2:]
        for arg in args:
            flag = 0
            for k in prof1.data.iterkeys():
                if arg == k.split()[0]:
                    flag = 1
            if not flag:
                print "ERROR: ", arg, " is not a valid parameter.  Type -options after -plot for help."
                sys.exit()
            prof1.plot(arg)
        plt.show()
    else:
        print "ERROR: Invalid number of arguments.  Type profiles_gen for help."
        sys.exit()
if len(sys.argv) <= 2:
    print "ERROR: Please specify data to be plotted, or ask for"
    print "options."
