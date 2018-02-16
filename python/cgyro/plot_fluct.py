import sys
import numpy as np
import os.path
from matplotlib import rc
import matplotlib.pyplot as plt
from gacodefuncs import *
from cgyro.data import cgyrodata
import gapy

# Use first 3 args to define plot and font size 
rc('text',usetex=True)
rc('font',size=int(sys.argv[12]))
#data_in.lx = int(sys.argv[2])
#data_in.ly = int(sys.argv[3])

ftype = sys.argv[1]
moment = sys.argv[2]
species = int(sys.argv[3])
ymin = sys.argv[4]
ymax = sys.argv[5]
nx = int(sys.argv[6])
ny = int(sys.argv[7])
istr = sys.argv[8]
fmin = sys.argv[9]
fmax = sys.argv[10]
colormap = sys.argv[11]

sim = cgyrodata('./')
nt = sim.n_time
nr = sim.n_radial
nn = sim.n_n
ns = sim.n_species

epx = np.zeros([nx,nr],dtype=np.complex)
eny = np.zeros([ny,nn],dtype=np.complex)
x = np.zeros([nx])
y = np.zeros([ny])

if istr == '-1':
    ivec = range(nt)
else:
    ivec = str2list(istr)

# Set filename and title
if (moment == 'n'):
    fdata = 'out.cgyro.kxky_n'
    title = r'$\delta \mathrm{n}$'
elif (moment == 'e'):
    fdata = 'out.cgyro.kxky_e'
    title = r'$\delta \mathrm{E}$'
elif (moment == 'phi'):
    fdata = 'out.cgyro.kxky_phi'
    title = r'$\delta\phi$'

# Check to see if it exists
if not os.path.isfile(fdata):
    print fdata+' does not exist.  Try -moment phi'
    sys.exit()

# WARNING: Assumes theta_plot=1 
#(2,n_radial,n_species,n_n,nt)
if (moment == 'phi'):
    n_chunk = 2*nr*nn
else:
    n_chunk = 2*nr*ns*nn
    
m = 0
i = 0

aa = np.zeros([n_chunk])
for line in open(fdata):
    aa[m] = float(line)
    m = m+1
    if m == n_chunk:
        i = i+1
        m = 0
        print 'INFO: (plot_fluct) Time index '+str(i) 
        if i in ivec:
            if (moment == 'phi'):
                a = np.reshape(aa,(2,nr,nn),order='F')
                c = a[0,:,:]+1j*a[1,:,:]
            else:
                a = np.reshape(aa,(2,nr,ns,nn),order='F')
                c = a[0,:,species,:]+1j*a[1,:,species,:]
                
            f = np.zeros([nx,ny],order='F')
            gapy.realfluct(c,f)
            if fmin == 'auto':
                f0=np.min(f)
                f1=np.max(f)
            else:
                f0=float(fmin)
                f1=float(fmax)

            xp = x/(2*np.pi)*sim.length
            yp = y/sim.ky[1]
            aspect = max(yp)/max(xp)

            fig = plt.figure(figsize=(10,10*aspect))
            ax = fig.add_subplot(111)
            ax.set_title(title)
            ax.set_xlabel(r'$x/\rho_s$')
            ax.set_ylabel(r'$y/\rho_s$')
            ax.set_aspect('equal')
        
            levels = np.arange(f0,f1,(f1-f0)/256)
            ax.contourf(xp,yp,np.transpose(f),levels,cmap=plt.get_cmap(colormap))
            print 'INFO: (plot_fluct) min,max = ',f0,f1

            fig.tight_layout(pad=0.5)
            if ftype == 'screen':
                plt.show()
            else:
                fname = fdata+str(i)
                # Filename uses frame number 
                plt.savefig(str(i)+'.png')
                # Close each time to prevent memory accumulation
                plt.close()

            if i == max(ivec):
                sys.exit()





