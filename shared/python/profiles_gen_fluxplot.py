import sys
import matplotlib.pyplot as plt
from profiles_genData import profiles_genData
from math import *

prof1 = profiles_genData(sys.argv[1])
m1 = float(sys.argv[3])
mteq = []
feq = []

fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.set_ylabel('Z (m)')
ax.set_xlabel('R (m)')
ax.set_title('Flux Surfaces')

if sys.argv[2] == '-m':
    if len(sys.argv) == 4:
        if m1 > 1:
            print "ERROR: Rho cannot be more than 1."
            sys.exit()
        if m1 < 0:
            print "ERROR: Rho cannot be less than 0."
            sys.exit()
        m1 = int(floor(m1 * 50))
        mteq.append(prof1.compute_mtypeeq(m1))
        ax.plot(mteq[0][0], mteq[0][1])
        plt.show()
    elif len(sys.argv) == 6:
        m2 = float(sys.argv[4])
        min1 = min(m1, m2)
        max1 = max(m1, m2)
        step = float(sys.argv[5])
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
            print "Changing step to 0.02"
        min1 = int(floor(min1 * 50))
        max1 = int(ceil(max1 * 50))
        step = int(floor(step * 50))
        for r in range(min1, max1 + 1, step):
            mteq.append(prof1.compute_mtypeeq(r))
            ax.plot(mteq[(r - min1) / step][0], mteq[(r - min1) / step][1], 'b')
        plt.show()
    else:
        print "Strange number of parameters.  Type profiles_gen for help."

if sys.argv[2] == '-f':
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
        plt.show()
    elif len(sys.argv) == 6:
        m2 = float(sys.argv[4])
        step = float(sys.argv[5])
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
        plt.show()
    else:
        print "Strange number of parameters.  Type profiles_gen for help."

if sys.argv[2] == '-c':
    if len(sys.argv) == 4:
        if m1 > 1:
            print "ERROR: Rho cannot be more than 1."
            sys.exit()
        if m1 < 0:
            print "ERROR: Rho cannot be less than 0." 
            sys.exit()
        m1 = int(floor(m1 * 50))
        feq.append(prof1.compute_fouriereq(m1))
        ax.plot(feq[0][0], feq[0][1])
        mteq.append(prof1.compute_mtypeeq(m1))
        ax.plot(mteq[0][0], mteq[0][1])
        plt.show()
    elif len(sys.argv) == 6:
        m2 = float(sys.argv[4])
        step = float(sys.argv[5])
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
        plt.show()
    else:
        print "Strange number of parameters.  Type profiles_gen for help."
