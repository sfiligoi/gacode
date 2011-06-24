"""profiles_gen_fluxplot.py contains plotting routines for data parsed with profiles_genData."""

import sys
import matplotlib.pyplot as plt
from profiles_genData import profiles_genData
from math import *

if len(sys.argv) < 4:
    print "ERROR: Too few arguments.  Type profiles_gen for help."
    sys.exit()
if len(sys.argv) > 6:
    print "ERROR: Too many arguments.  Type profiles_gen for help."
try:
    prof1 = profiles_genData(sys.argv[1])
except IOError:
    print sys.argv[1]
    print "does not contain file input.profiles and/or file input.profiles.geo."
    print "Type profiles_gen for help."
    sys.exit()

mteq = []
feq = []

#Produces a matplotlib figure object and creates the labels.
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.set_ylabel('Z (m)')
ax.set_xlabel('R (m)')

#Checks to see which type of plot is desired: Miller-type, Fourier, or a
#comparison of the two.
if sys.argv[2] == '-m':
    try:
        m1 = float(sys.argv[3])
    except ValueError:
        print sys.argv[3] + " is not a vaild number."
        sys.exit()
#   Checks to see if only one radius is desired to be plotted, or more.
    if len(sys.argv) == 4:
#       Bounds checking
        if m1 > 1:
            print "ERROR: Rho cannot be more than 1."
            sys.exit()
        if m1 < 0:
            print "ERROR: Rho cannot be less than 0."
            sys.exit()
#       Converts m1 to proper index number for getting data out of prof1.
        m1 = int(floor(m1 * 50))
#       Computes fits and stores them in mteq.
        mteq.append(prof1.compute_mtypeeq(m1))
        ax.plot(mteq[0][0], mteq[0][1], 'b')
        ax.set_xlim(0.55, 2.75)
        ax.set_ylim(-1.0, 1.2)
        plt.show()
    elif len(sys.argv) == 6:
        try:
            m2 = float(sys.argv[4])
        except ValueError:
            print sys.argv[4] + " is not a vaild number."
            sys.exit()
        try:
            step = float(sys.argv[5])
        except ValueError:
            print sys.argv[5] + " is not a vaild number."
            sys.exit()
        min1 = min(m1, m2)
        max1 = max(m1, m2)
        if min1 >= 1:
            print "ERROR: Min must be less than 1."
            sys.exit()
        if min1 < 0:
            print "ERROR: Min cannot be less than 0."
            sys.exit()
        if max1 > 1:
            print "ERROR: Max cannot be greater than 1."
            sys.exit()
        if max1 <= 0:
            print "ERROR: Max must be greater than 0."
            sys.exit()
        if step < prof1.step:
            step = prof1.step
            print "WARNING: " + str(prof1.step) + " is smallest step size available."
            print "Changing step size to " + str(prof1.step) + "."
        min1 = int(floor(min1 * 50))
        max1 = int(ceil(max1 * 50))
        step = int(floor(step * 50))
        for r in range(min1, max1 + 1, step):
            mteq.append(prof1.compute_mtypeeq(r))
            ax.plot(mteq[(r - min1) / step][0], mteq[(r - min1) / step][1], 'b')
        ax.set_xlim(0.55, 2.75)
        ax.set_ylim(-1.0, 1.2)
        ax.set_title('Flux Surfaces with Min = ' + str(min1/50.0) + ', Max = ' + str(max1/50.0) + ', and Step Size = ' + str(step/50.0))
        plt.show()
    else:
        print "Strange number of parameters.  Type profiles_gen for help."
elif sys.argv[2] == '-f':
    try:
        m1 = float(sys.argv[3])
    except ValueError:
        print sys.argv[3] + " is not a vaild number."
        sys.exit()
    if len(sys.argv) == 4:
        if m1 > 1:
            print "ERROR: Rho cannot be more than 1."
            sys.exit()
        if m1 < 0:
            print "ERROR: Rho cannot be less than 0." 
            sys.exit()
        m1 = int(floor(m1 * 50))
        feq.append(prof1.compute_fouriereq(m1))
        ax.plot(feq[0][0], feq[0][1], 'r')
        ax.set_xlim(0.55, 2.75)
        ax.set_ylim(-1.0, 1.2)
        plt.show()
    elif len(sys.argv) == 6:
        try:
            m2 = float(sys.argv[4])
        except ValueError:
            print sys.argv[4] + " is not a vaild number."
            sys.exit()
        try:
            step = float(sys.argv[5])
        except ValueError:
            print sys.argv[5] + " is not a vaild number."
            sys.exit()
        min1 = min(m1, m2)
        max1 = max(m1, m2)
        if min1 >= 1:
            print "ERROR: Min must be less than 1."
            sys.exit()
        if min1 < 0:
            print "ERROR: Min cannot be less than 0."
            sys.exit()
        if max1 > 1:
            print "ERROR: Max cannot be greater than 1."
            sys.exit()
        if max1 <= 0:
            print "ERROR: Max must be greater than 0."
            sys.exit()
        if step < 0.02:
            step = 0.02
            print "WARNING: 0.02 is smallest step size available."
            print "Changing step size to 0.02"
        min1 = int(floor(min1 * 50))
        max1 = int(ceil(max1 * 50))
        step = int(floor(step * 50))
        for r in range(min1, max1 + 1, step):
            feq.append(prof1.compute_fouriereq(r))
            ax.plot(feq[(r - min1) / step][0], feq[(r - min1) / step][1], 'r')
        ax.set_xlim(0.55, 2.75)
        ax.set_ylim(-1.0, 1.2)
        plt.show()
    else:
        print "Strange number of parameters.  Type profiles_gen for help."
elif sys.argv[2] == '-c':
    try:
        m1 = float(sys.argv[3])
    except ValueError:
        print sys.argv[3] + " is not a vaild number."
    if len(sys.argv) == 4:
        if m1 > 1:
            print "ERROR: Rho cannot be more than 1."
            sys.exit()
        if m1 < 0:
            print "ERROR: Rho cannot be less than 0." 
            sys.exit()
        m1 = int(floor(m1 * 50))
        feq.append(prof1.compute_fouriereq(m1))
        ax.plot(feq[0][0], feq[0][1], 'r')
        mteq.append(prof1.compute_mtypeeq(m1))
        ax.plot(mteq[0][0], mteq[0][1], 'b')
        ax.set_xlim(0.55, 2.75)
        ax.set_ylim(-1.0, 1.2)
        plt.show()
    elif len(sys.argv) == 6:
        try:
            m2 = float(sys.argv[4])
        except ValueError:
            print sys.argv[4] + " is not a vaild number."
        try:
            step = float(sys.argv[5])
        except ValueError:
            print sys.argv[5] + " is not a vaild number."
        min1 = min(m1, m2)
        max1 = max(m1, m2)
        if min1 >= 1:
            print "ERROR: Min must be less than 1."
            sys.exit()
        if min1 < 0:
            print "ERROR: Min cannot be less than 0."
            sys.exit()
        if max1 > 1:
            print "ERROR: Max cannot be greater than 1."
            sys.exit()
        if max1 <= 0:
            print "ERROR: Max must be greater than 0."
            sys.exit()
        if step < 0.02:
            step = 0.02
            print "WARNING: 0.02 is smallest step size available."
            print "Changing step size to 0.02"
        min1 = int(floor(min1 * 50))
        max1 = int(ceil(max1 * 50))
        step = int(floor(step * 50))
        for r in range(min1, max1 + 1, step):
            feq.append(prof1.compute_fouriereq(r))
            ax.plot(feq[(r - min1) / step][0], feq[(r - min1) / step][1], 'r')
            mteq.append(prof1.compute_mtypeeq(r))
            ax.plot(mteq[(r - min1) / step][0], mteq[(r - min1) / step][1], 'b')
            ax.set_xlim(0.55, 2.75)
            ax.set_ylim(-1.0, 1.2)
        plt.show()
    else:
        print "Strange number of parameters.  Type profiles_gen for help."
else:
    print "ERROR: Incorrect plot type.  Type profiles_gen for help."
