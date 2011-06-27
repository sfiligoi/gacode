"""profiles_gen_fluxplot.py contains plotting routines for data parsed with profiles_genData."""

import sys
import matplotlib.pyplot as plt
from profiles_genData import profiles_genData
from math import *

if len(sys.argv) < 4:
    print "ERROR: Too few arguments.  Type profiles_gen for help."
    sys.exit()

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
        print "ERROR: " + sys.argv[3] + " is not a vaild number."
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

#       Computes fits and stores them in mteq.
        mteq.append(prof1.compute_mtypeeq(m1))
        r = mteq[0][0]
        z = mteq[0][1]
        rmaj = float(mteq[0][2])
        zmag = float(mteq[0][3])
        ax.plot(r, z, 'b')
        aspect = max((max(z) - min(z) + .5), (max(r) - min(r) + .5))
        ax.set_xlim(rmaj - aspect/2, rmaj + aspect/2)
        ax.set_ylim(zmag - aspect/2, zmag + aspect/2)
        ax.set_title('Flux Surface at Rho= ' + str(prof1.get('rho')[prof1.match(m1, prof1.get('rho'))]))
        plt.show()

    elif len(sys.argv) == 6:
        try:
            m2 = float(sys.argv[4])
        except ValueError:
            print "ERROR: " + sys.argv[4] + " is not a vaild number."
            sys.exit()
        try:
            n = float(sys.argv[5])
        except ValueError:
            print "ERROR: " + sys.argv[5] + " is not a vaild number."
            sys.exit()
        if n != floor(n):
            print "ERROR: Number of flux surfaces must be an integer."
            sys.exit()
        if n < 2:
            print "ERROR: Number of flux surfaces must at least 2."
            sys.exit()

        min1 = min(m1, m2)
        inc = min1
        max1 = max(m1, m2)
        step = (max1 - min1) / n
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

        while inc < max1:
            mteq.append(prof1.compute_mtypeeq(inc))
            r = mteq[0][0]
            z = mteq[0][1]
            rmaj = float(mteq[0][2])
            zmag = float(mteq.pop()[3])
            ax.plot(r, z, 'b')
            inc = inc + step
        aspect = max((max(z) - min(z) + .5), (max(r) - min(r) + .5))
        ax.set_xlim(rmaj - aspect/2, rmaj + aspect/2)
        ax.set_ylim(zmag - aspect/2, zmag + aspect/2)
        ax.set_title(str(int(n)) + ' Flux Surfaces between ' + str(min1) + ' and ' + str(max1) + '.')
        plt.show()

    else:
        print "Strange number of parameters.  Type profiles_gen for help."



elif sys.argv[2] == '-f':
    try:
        m1 = float(sys.argv[3])
    except ValueError:
        print "ERROR: " + sys.argv[3] + " is not a vaild number."
        sys.exit()

    if len(sys.argv) == 4:
        if m1 > 1:
            print "ERROR: Rho cannot be more than 1."
            sys.exit()
        if m1 < 0:
            print "ERROR: Rho cannot be less than 0." 
            sys.exit()

        feq.append(prof1.compute_fouriereq(m1))
        r = feq[0][0]
        z = feq[0][1]
        rmaj = float(feq[0][2])
        zmag = float(feq[0][3])
        ax.plot(r, z, 'r')
        aspect = max((max(z) - min(z) + .5), (max(r) - min(r) + .5))
        ax.set_xlim(rmaj - aspect/2, rmaj + aspect/2)
        ax.set_ylim(zmag - aspect/2, zmag + aspect/2)
        ax.set_title('Flux Surface at Rho= ' + str(prof1.get('rho')[prof1.match(m1, prof1.get('rho'))]))
        plt.show()

    elif len(sys.argv) == 6:
        try:
            m2 = float(sys.argv[4])
        except ValueError:
            print "ERROR: " + sys.argv[4] + " is not a vaild number."
            sys.exit()
        try:
            n = float(sys.argv[5])
        except ValueError:
            print "ERROR: " + sys.argv[5] + " is not a vaild number."
            sys.exit()
        if n != floor(n):
            print "ERROR: Number of flux surfaces must be an integer."
            sys.exit()
        if n < 2:
            print "ERROR: Number of flux surfaces must at least 2."
            sys.exit()

        min1 = min(m1, m2)
        inc = min1
        max1 = max(m1, m2)
        step = (max1 - min1) / n
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

        while inc < max1:
            feq.append(prof1.compute_fouriereq(inc))
            r = feq[0][0]
            z = feq[0][1]
            rmaj = float(feq[0][2])
            zmag = float(feq.pop()[3])
            ax.plot(r, z, 'r')
            inc = inc + step
        aspect = max((max(z) - min(z) + .5), (max(r) - min(r) + .5))
        ax.set_xlim(rmaj - aspect/2, rmaj + aspect/2)
        ax.set_ylim(zmag - aspect/2, zmag + aspect/2)
        ax.set_title(str(int(n)) + ' Flux Surfaces between ' + str(min1) + ' and ' + str(max1) + '.')
        plt.show()

    else:
        print "Strange number of parameters.  Type profiles_gen for help."



elif sys.argv[2] == '-c':
    try:
        m1 = float(sys.argv[3])
    except ValueError:
        print "ERROR: " + sys.argv[3] + " is not a vaild number."
        sys.exit()

    if len(sys.argv) == 4:
        if m1 > 1:
            print "ERROR: Rho cannot be more than 1."
            sys.exit()
        if m1 < 0:
            print "ERROR: Rho cannot be less than 0." 
            sys.exit()

        feq.append(prof1.compute_fouriereq(m1))
        fr = feq[0][0]
        fz = feq[0][1]
        ax.plot(fr, fz, 'r')
        mteq.append(prof1.compute_mtypeeq(m1))
        mr = mteq[0][0]
        mz = mteq[0][1]
        rmaj = float(mteq[0][2])
        zmag = float(mteq[0][3])
        ax.plot(mr, mz, 'b')
        aspect = max((max(fz) - min(fz) + .5), (max(fr) - min(fr) + .5), (max(mz) - min(mz) + .5), (max(mr) - min(mr) + .5))
        ax.set_xlim(rmaj - aspect/2, rmaj + aspect/2)
        ax.set_ylim(zmag - aspect/2, zmag + aspect/2)
        ax.set_title('Flux Surface at Rho= ' + str(prof1.get('rho')[prof1.match(m1, prof1.get('rho'))]))
        plt.show()

    elif len(sys.argv) == 6:
        try:
            m2 = float(sys.argv[4])
        except ValueError:
            print "ERROR: " + sys.argv[4] + " is not a vaild number."
            sys.exit()
        try:
            n = float(sys.argv[5])
        except ValueError:
            print "ERROR: " + sys.argv[5] + " is not a vaild number."
            sys.exit()
        if n != floor(n):
            print "ERROR: Number of flux surfaces must be an integer."
            sys.exit()
        if n < 2:
            print "ERROR: Number of flux surfaces must at least 2."
            sys.exit()

        min1 = min(m1, m2)
        inc = min1
        max1 = max(m1, m2)
        step = (max1 - min1) / n
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

        while inc < max1:
            feq.append(prof1.compute_fouriereq(inc))
            fr = feq[0][0]
            fz = feq.pop()[1]
            ax.plot(fr, fz, 'r')
            mteq.append(prof1.compute_mtypeeq(inc))
            mr = mteq[0][0]
            mz = mteq[0][1]
            rmaj = float(mteq[0][2])
            zmag = float(mteq.pop()[3])
            ax.plot(mr, mz, 'b')
            inc = inc + step
        aspect = max((max(fz) - min(fz)), (max(fr) - min(fr)), (max(mz) - min(mz)), (max(mr) - min(mr)))
        ax.set_xlim(rmaj - aspect/2, rmaj + aspect/2)
        ax.set_ylim(zmag - aspect/2, zmag + aspect/2)
        ax.set_title(str(int(n)) + ' Flux Surfaces between ' + str(min1) + ' and ' + str(max1) + '.')
        plt.show()

    else:
        print "Strange number of parameters.  Type profiles_gen for help."


else:
    print "ERROR: Incorrect plot type.  Type profiles_gen for help."
