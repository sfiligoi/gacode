import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
w = float(sys.argv[2])
ymax  = sys.argv[3]
nstr  = sys.argv[4]

sim = cgyrodata('./')
sim.getbigfield()

fig = plt.figure(figsize=(12,6))

color = ['m','k','b','c']

kx  = sim.kx
ave = np.zeros(sim.n_radial)

imin=iwindow(sim.t,w)

xlabel=r'$k_x \rho_s$'
    
dk = kx[1]-kx[0]
x0 = kx[-1]+dk

ax = fig.add_subplot(1,1,1)
ax.set_xlabel(xlabel)
ax.set_ylabel(r'$\phi_\mathrm{RMS}$',color='k')
ax.set_title(r'$'+str(sim.t[imin])+' < (c_s/a) t < '+str(sim.t[-1])+'$')

if nstr == 'null':
    y = np.sum(sim.phisq,axis=1)
    for j in range(sim.n_radial):
        ave[j] = average(y[j,:],sim.t,w)
    ax.plot(kx,np.sqrt(ave[:]),color=color[0],alpha=0.5)
else:
    y = np.zeros([sim.n_radial,sim.n_time])
    nvec = str2list(nstr)
    print 'INFO: (plot_kx_phi) n = '+str(nvec)
    for n in nvec:
        num = r'$n='+str(n)+'$'
        y[:] = sim.phisq[:,n,:]
        for j in range(sim.n_radial):
            ave[j] = average(sim.phisq[j,n,:],sim.t,w)
        ax.plot(kx+dk/2,np.sqrt(ave[:]),ls='steps',label=num)
            
ax.set_xlim([-x0,x0])
ax.set_yscale('log')
ax.set_ylim([1e-5,1.0])
             
ax.plot(kx,sim.radialdiss*ax.get_ylim()[1]*0.5,linewidth=2,color='k',alpha=0.2)


if sim.n_n > 16:
    ax.legend(loc=4, ncol=5, prop={'size':12})
else:
    ax.legend(loc=4, ncol=6, prop={'size':12})

if ftype == 'screen':
   plt.show()
else:
    fname=fname+'.'+ftype
    print 'INFO: (plot_kx_phi) Created '+fname
    plt.savefig(fname)
