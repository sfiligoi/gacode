import sys
import numpy as np
import matplotlib.pyplot as plt
from pyrats.data import GYROData

directory=sys.argv[1]
n_field=sys.argv[3]
try:
    n_kinetic=int(sys.argv[2])
    n_moment=int(sys.argv[4])
except ValueError:
    print "ERROR: Invalid parameter.  Try again."
    sys.exit()

try:
    sim1 = GYROData(directory)
except IOError:
    print "ERROR: Invalid directory location.  Try again."
    sys.exit()
if (n_kinetic > sim1.profile['n_kinetic']) or (n_moment > 2):
    print "ERROR: Parameter out of range.  Try again."
    sys.exit()
sim1.make_gbflux()
fig = plt.figure(1)
ax = fig.add_subplot(111)
if n_field == 's':
    ax.plot(np.sum(sim1.gbflux, axis=1)[n_kinetic, n_moment, :])
    ax.set_title("Species " + str(n_kinetic) + ", field " + str(n_field) + ", moment " + str(n_moment))
else:
     try:
        n_field = int(n_field)
        if (n_field > sim1.profile['n_field']):
            print "ERROR: Parameter out of range.  Try again."
            sys.exit()
     except ValueError:
        print "ERROR: Invalid parameter.  Try again."
        sys.exit()
     ax.plot(sim1.gbflux[n_kinetic, n_field, n_moment, :])
     ax.set_title("Species " + str(n_kinetic) + ", field " + str(n_field) + ", moment " + str(n_moment))
ax.set_xlabel("Time")
if n_moment == 0:
    ax.set_ylabel("GAMMA_" + str(n_kinetic) + "/GAMMA_{GB}")
elif n_moment == 1:
    ax.set_ylabel("Q_{0," + str(n_kinetic) + "}/Q_{GB}")
elif n_moment == 2:
    ax.set_ylabel("PI_" + str(n_kinetic) + "/PI_{GB}")
elif n_moment == 3:
    ax.set_ylabel("S^tur_{W," + str(n_kinetic) + "}/S_{GB}")
plt.show()
