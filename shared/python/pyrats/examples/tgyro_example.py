from pyrats.tgyro.data import TGYROData
from pyrats.neo.data import NEOData
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager

prop = matplotlib.font_manager.FontProperties(size=16)

path   = '/home/buuck/'
tgyro1 = TGYROData(path+'treg01')
neo1   = NEOData(path+'treg01')

rtext = r'$r/a$'

fig = plt.figure(figsize=(14,10))
fig.subplots_adjust(wspace=0.3)
fig.subplots_adjust(left=0.07,right=0.93)
fig.subplots_adjust(top=0.9,bottom=0.12)

#-------------------------------------------------
# Construct primitive NEO data from TGLF* dirs.
n=12
qi = np.zeros(n+1)
qe = np.zeros(n+1)
te = neo1.control['T0_over_Tnorm'].data[:, 1]
rhostar = neo1.control['rho_star'].data
qgb = te**2.5*rhostar**2
qi[1:] = neo1.transport['Q'].data[:, 0]/qgb
qe[1:] = neo1.transport['Q'].data[:, 1]/qgb
#-------------------------------------------------

# Turbulence
ax = fig.add_subplot(221)
ax.set_ylabel(r'Turbulent Fluxes (GB units)')

ax.set_xlabel(rtext,fontsize=16)

a=ax.plot(tgyro1.data.get('r/a')[-1], tgyro1.data.get('eflux_i_tur')[-1],'b')
b=ax.plot(tgyro1.data.get('r/a')[-1], tgyro1.data.get('eflux_e_tur')[-1],'r')
ax.legend((a,b),(r'$Q_D/Q_{\rm GB}$',
                 r'$Q_e/Q_{\rm GB}$'),'upper left',prop=prop)

# Neoclassical
ax = fig.add_subplot(222)

ax.set_xlabel(rtext,size=16)
ax.set_ylabel('Neoclassical Fluxes (GB units)')

print sorted(tgyro1.data.keys())
a=ax.plot(tgyro1.data.get('r/a')[-1], tgyro1.data.get('eflux_i_neo')[-1],'b')
b=ax.plot(tgyro1.data.get('r/a')[-1], tgyro1.data.get('eflux_e_neo')[-1],'r')
c=ax.plot(tgyro1.data.get('r/a')[-1], qi,'k--')
d=ax.plot(tgyro1.data.get('r/a')[-1], qe,'k--')
ax.legend((a,b,c,d),(r'$Q_D/Q_{\rm GB}$',
                     r'$Q_e/Q_{\rm GB}$',
                     r'$Q_D/Q_{\rm GB}$',
                     r'$Q_e/Q_{\rm GB}$'),'upper left',prop=prop)

#================================================================

qgb = tgyro1.data.get('Q_GB')[4]

# Turbulence
ax = fig.add_subplot(223)

ax.set_xlabel(rtext,size=16)
ax.set_ylabel(r'Turbulent Fluxes (MW/m$^2$)')

a=ax.plot(tgyro1.data.get('r/a')[-1], tgyro1.data.get('eflux_i_tur')[-1]*qgb,'b')
b=ax.plot(tgyro1.data.get('r/a')[-1], tgyro1.data.get('eflux_e_tur')[-1]*qgb,'r')
ax.legend((a,b),(r'$Q_D/Q_{\rm GB}$',
                 r'$Q_e/Q_{\rm GB}$'),'upper left',prop=prop)

# Neoclassical
ax = fig.add_subplot(224)

ax.set_xlabel(rtext,size=16)
ax.set_ylabel(r'Neoclassical Fluxes (MW/m$^2$)')

a=ax.plot(tgyro1.data.get('r/a')[-1], tgyro1.data.get('eflux_i_neo')[-1]*qgb,'b')
b=ax.plot(tgyro1.data.get('r/a')[-1], tgyro1.data.get('eflux_e_neo')[-1]*qgb,'r')
ax.legend((a,b),(r'$Q_D/Q_{\rm GB}$',
                 r'$Q_e/Q_{\rm GB}$'),'upper left',prop=prop)

print tgyro1.data.keys()

plt.show()
