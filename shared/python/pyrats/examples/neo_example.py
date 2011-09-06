from pyrats.neo.data import NEOData
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager
prop = matplotlib.font_manager.FontProperties(size=16)

path = '/home/candy/mysim/pyrats/'
neo1 = NEOData(path+'neo1')

rtext = r'$r\,({\rm m})$'

# NOTE:
#
# This is an radial scan, so
#
# m_norm = mD
# n_norm = n_{first species} = nD
# T_norm = T_{first species} = TD

fig = plt.figure(figsize=(14,5))
fig.subplots_adjust(wspace=0.3)
fig.subplots_adjust(left=0.07,right=0.93)
fig.subplots_adjust(top=0.9,bottom=0.12)

#-----------------------------------------------------------
ax = fig.add_subplot(131)

ax.set_xlabel(rtext,size=16)
ax.set_ylabel(r'$Q_{\rm GB}$',size=16)

# For this case: 0=D, 1=C, 2=e:
ne = neo1.equil['n'][:,2]
te = neo1.equil['T'][:,2]
rhostar = neo1.equil['rho_star']

qgb = ne*te**2.5*rhostar**2

q = neo1.transport['Q']
r = neo1.transport_exp['r']

# (0=D,1=e,2=C)
a=ax.plot(r,q[:,0]/qgb,'b')
b=ax.plot(r,q[:,1]/qgb,'r')
c=ax.plot(r,q[:,2]/qgb,'k')

ax.legend((a,b,c),(r'$Q_D$',
                   r'$Q_C$',
                   r'$Q_e$'),'upper left',prop=prop)
#-----------------------------------------------------------

#-----------------------------------------------------------
ax = fig.add_subplot(132)

ax.set_xlabel(rtext,size=16)
ax.set_ylabel(r'${\rm km/s}$',size=16)

v0   = neo1.equil['v_norm_over_a']/1000
vtor = neo1.transport['vphi']

# (0=D,1=e,2=c)
a=ax.plot(r,vtor[:,0]*v0,'b')
b=ax.plot(r,vtor[:,1]*v0,'r')
c=ax.plot(r,vtor[:,2]*v0,'k')

ax.legend((a,b,c),(r'$v_D$',
                   r'$v_C$',
                   r'$v_e$'),'upper left',prop=prop)
#-----------------------------------------------------------

#-----------------------------------------------------------
ax = fig.add_subplot(133)

ax.set_xlabel(rtext,size=16)
ax.set_ylabel(r'$\langle j_\parallel B \rangle/(e n_e v_D B_{\rm unit})$',size=16)

nd = neo1.equil['n'][:,0]
jb = neo1.transport['jparB']*nd/ne

a=ax.plot(r,jb,'b')
#-----------------------------------------------------------

plt.show()
