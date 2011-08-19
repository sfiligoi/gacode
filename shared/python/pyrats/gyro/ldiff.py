"""This file is executed by the bash script gyro_plot when a listing of the
diffusivities is requested."""

from pyrats.gyro.data import GYROData
import sys
import numpy as np

sim       = GYROData(sys.argv[1])
field     = sys.argv[2]
i_moment  = int(sys.argv[3])

n_field   = int(sim.profile['n_field'])
n_kinetic = int(sim.profile['n_kinetic'])

sim.make_diff()

t    = sim.t['(cbar_s/a)t']
flux = sim.diff

# b is collection of all arrays to be plotted
b = np.zeros((len(t),n_kinetic+1))

b[:,0] = t

# Manage field
if field == 's':
    flux0 = np.sum(flux,axis=1)
    ftag = 'TOT   '
else:
    i_field = int(field)
    flux0 = flux[:,i_field,:,:]
    if i_field == 0: 
        ftag = 'ES    '
    if i_field == 1: 
        ftag = 'EM    '
    if i_field == 2: 
        ftag = 'COM   '

# Manage moment
if i_moment == 0: 
    mtag = 'D [GB]      '
if i_moment == 1: 
    mtag = 'CHI [GB]    '

ul = '----------  '

line1 = '                '
line2 = '    (cs/a)t     '
line3 = '    '+ul

# Manage species
for i in range(n_kinetic):
    b[:,i+1] = flux0[i,i_moment,:]
    if i == n_kinetic-1:
        stag = 'elec  '
    else:
        stag = 'ion-'+str(i)+' '

    line1 = line1+mtag
    line2 = line2+stag+ftag
    line3 = line3+ul

np.set_printoptions(precision=3,suppress=False,threshold=100000)

print line1
print line2
print line3
print b
