import sys
import numpy as np
from gacodeplotdefs import *
from gyro.data import GYROData

sim       = GYROData(sys.argv[1])
field     = sys.argv[2]
i_moment  = int(sys.argv[3])
window    = float(sys.argv[4])
ftype     = sys.argv[5]
datafile  = sys.argv[6]

n_field   = int(sim.profile['n_field'])
n_kinetic = int(sim.profile['n_kinetic'])
n_n       = int(sim.profile['n_n'])

# Need to read gbflux_n data
sim.read_gbflux_n()

t    = sim.t['(c_s/a)t']
flux = sim.gbflux_n

# Manage field
if field == 's':
    flux0 = np.sum(flux,axis=1)
    ftag  = sim.tagfield[3]
else:
    i_field = int(field)
    flux0 = flux[:,i_field,:,:,:]
    ftag  = sim.tagfield[i_field]

# Manage moment
mtag = sim.tagmom[i_moment]

#======================================
fig = plt.figure(figsize=(7*n_kinetic,6))
#=====================================

color = ['k','m','b','c']

k   = sim.profile['kt_rho']
dk  = k[1]-k[0]
ave = np.zeros((n_n))

# Determine tmin
imin=0
for i in range(len(t)):
    if t[i] < (1.0-window)*t[len(t)-1]:
        imin = i+1

if datafile == 'none':

    # Plot data to screen or image file.

    for i in range(n_kinetic):
        ax = fig.add_subplot(1,n_kinetic,i+1)
        ax.set_xlabel(r'$k_\theta \rho_s$',fontsize=GFONTSIZE)
        ax.set_ylabel(r'$'+mtag+' \;('+ftag+')$',color='k',fontsize=GFONTSIZE)
        stag = sim.tagspec[i]
        for j in range(n_n):
            ave[j] = average(flux0[i,i_moment,j,:],t,window)
        ax.set_title(stag+r': $'+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+'$',fontsize=GFONTSIZE)
        ax.bar(k-dk/2.0,ave,width=dk/1.1,color=color[i],alpha=0.4,edgecolor='black')
        
    if ftype == 'screen':
        plt.show()
    else:
        outfile = 'gbflux_n.'+ftype
        plt.savefig(outfile)

else:

    # Write data to datafile

    arr = np.zeros([len(k),n_kinetic+1])
    arr[:,0] = k
    stag = '# (k_y rho_s'
    for i in range(n_kinetic):
        for j in range(n_n):
            ave[j] = average(flux0[i,i_moment,j,:],t,window)
        arr[:,i+1] = ave
        stag = stag+' , '+sim.tagspec[i]
            
    fid = open(datafile,'w')
    fid.write('# Moment  : '+mtag+'\n')
    fid.write('# Field   : '+ftag+'\n')
    fid.write('# Time    : '+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+'\n')
    fid.write(stag+')\n')
    np.savetxt(fid,arr,fmt='%.5e')
    fid.close()

