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
from profiles_genData import profiles_genData
from profiles_generrorcheck import *
try:
    import matplotlib.pyplot as plt
except ImportError:
    print "This command requires matplotlib.  Please install matplotlib:"
    print "http://matplotlib.sourceforge.net/"
    sys.exit()
import math

args = ['ne', 'ni_1', 'Te', 'Ti_1']
n1 = 1
n2 = 1

#Error catching: no arguments mean that the default directory is used.
if len(sys.argv) == 1:
    prof1 = opendir('.')
    for arg in args:
        checkarg(prof1, arg)
        prof1.plot(arg, 2, 2)
    plt.show()

#If '-c' is the second argument, then this portion of code tries to open up the
#next two arguments as two different directories.
else:
    if sys.argv[1] == '-c':
        #And of course runs error checks when trying to open them
        prof1 = opendir(sys.argv[2])
        prof2 = opendir(sys.argv[3])
        #If there are more than 4 args, then the user is asking for more than
        #just the default behavior.
        if len(sys.argv) > 4:
            #If the next argument is '-options', then the program looks in both
            #directories and prints only those parameters which are nonzero in
            #both
            if sys.argv[4] == '-options':
                keys = []
                for k, v in prof1.data.iteritems():
                    s = 0
                    for item in v:
                        s = s + float(item)
                    if s != 0:
                        s = 0
                        for item in prof2.get(k):
                            s = s + float(item)
                        if s != 0:
                            keys.append(k.split()[0])
                keys.sort()
                for k in keys:
                    print k
            #If the next argument is a number, then the user must be trying to
            #set the dimensions
            elif sys.argv[4].isdigit():
                if len(sys.argv) > 6:
                    #So everything after must be parameters to plot
                    args = sys.argv[6:]
                n1 = setdim(sys.argv[4])
                n2 = setdim(sys.argv[5])
                #Check to make sure that the requested plots actually exist
                for arg in args:
                    checkarg(prof1, arg)
                    checkarg(prof2, arg)
                    #And finally we create the plots
                    prof1.plot(arg, n1, n2, 'orange')
                    prof2.plot(arg, n1, n2, 'b')
                plt.show()
        #If the 5th argument is not a digit, then it must be something to plot.
            elif sys.argv[4].isdigit() == False:
                args = sys.argv[4:]
                for arg in args:
                    checkarg(prof1, arg)
                    checkarg(prof2, arg)
                    prof1.plot(arg, n1, n2, 'orange')
                    prof2.plot(arg, n1, n2, 'b')
                plt.show()
            else:
                print "ERROR: Invalid number of arguments.  Type profiles_gen for help."
                sys.exit()
        #The default behavior is to plot whatever is in 'args' on a 2x2 plot.
        else:
            for arg in args:
                checkarg(prof1, arg)
                checkarg(prof2, arg)
                prof1.plot(arg, 2, 2, 'orange')
                prof2.plot(arg, 2, 2, 'b')
            plt.show()
#If the second argument is anything besides '-c', then it must be a directory location.
    else:
        #Unless its '-options'
        if sys.argv[1] == '-options':
            prof1 = opendir('.')
            getoptions(prof1)
        else:
            prof1 = opendir(sys.argv[1])
            if len(sys.argv) > 2:
                #This code segment prints all of the variables with nonzero
                #entries when '-options' is called.
                if sys.argv[2] == '-options':
                    getoptions(prof1)
                #More error catching; the rest of this is pretty similar to the
                #previous part.
                elif sys.argv[2].isdigit():
                    args = sys.argv[4:]
                    n1 = setdim(sys.argv[2])
                    n2 = setdim(sys.argv[3])
                    for arg in args:
                        checkarg(prof1, arg)
                        prof1.plot(arg, n1, n2)
                    plt.show()
                elif sys.argv[2].isdigit() == False:
                    args = sys.argv[2:]
                    for arg in args:
                        checkarg(prof1, arg)
                        prof1.plot(arg)
                    plt.show()
            else:
                for arg in args:
                    checkarg(prof1, arg)
                    prof1.plot(arg, 2, 2)
                plt.show()
