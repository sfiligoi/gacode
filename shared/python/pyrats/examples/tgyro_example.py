from pyrats.tgyro.data import TGYROData
from pyrats.neo.data import NEOData
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager

prop = matplotlib.font_manager.FontProperties(size=16)

path   = '/home/candy/mysim/pyrats/'
tgyro1 = TGYROData(path+'tgyro1')
neo1   = NEOData(path+'tgyro1')

rtext = r'$r/a$'

fig = plt.figure(figsize=(14,10))
fig.subplots_adjust(wspace=0.3)
fig.subplots_adjust(left=0.07,right=0.93)
fig.subplots_adjust(top=0.9,bottom=0.12)

#-------------------------------------------------
# Construct primitive NEO data from TGLF* dirs.
n = 12
qi = np.zeros(n+1)
qe = np.zeros(n+1)
for i in range(n):
    ppath = path+'tgyro1/TGLF'+str(i+1)
    te = neo1.control['T0_over_T_norm'].data[1]
    rhostar = neo1.control['rho_star'].data
    qgb = te**2.5*rhostar**2
    qi[i+1] = neo1.transport['Q'].data[0]/qgb
    qe[i+1] = neo1.transport['Q'].data[1]/qgb
#-------------------------------------------------

# Turbulence
ax = fig.add_subplot(221)
ax.set_ylabel(r'Turbulent Fluxes (GB units)')

ax.set_xlabel(rtext,fontsize=16)

a=ax.plot(tgyro1.get_r(), tgyro1.get_flux_i_turb(),'b')
b=ax.plot(tgyro1.get_r(), tgyro1.get_flux_e_turb(),'r')
ax.legend((a,b),(r'$Q_D/Q_{\rm GB}$',
                 r'$Q_e/Q_{\rm GB}$'),'upper left',prop=prop)

# Neoclassical
ax = fig.add_subplot(222)

ax.set_xlabel(rtext,size=16)
ax.set_ylabel('Neoclassical Fluxes (GB units)')

print tgyro1.flux_e.keys()
a=ax.plot(tgyro1.get_r(), tgyro1.get_flux_i_neo(),'b')
b=ax.plot(tgyro1.get_r(), tgyro1.flux_e['eflux_e_neo'][-1],'r')
c=ax.plot(tgyro1.get_r(), qi,'k--')
d=ax.plot(tgyro1.get_r(), qe,'k--')
ax.legend((a,b,c,d),(r'$Q_D/Q_{\rm GB}$',
                     r'$Q_e/Q_{\rm GB}$',
                     r'$Q_D/Q_{\rm GB}$',
                     r'$Q_e/Q_{\rm GB}$'),'upper left',prop=prop)

#================================================================

qgb = tgyro1.gyro_bohm_unit['Q_GB'][4]

# Turbulence
ax = fig.add_subplot(223)

ax.set_xlabel(rtext,size=16)
ax.set_ylabel(r'Turbulent Fluxes (MW/m$^2$)')

a=ax.plot(tgyro1.get_r(), tgyro1.get_flux_i_turb()*qgb,'b')
b=ax.plot(tgyro1.get_r(), tgyro1.get_flux_e_turb()*qgb,'r')
ax.legend((a,b),(r'$Q_D/Q_{\rm GB}$',
                 r'$Q_e/Q_{\rm GB}$'),'upper left',prop=prop)

# Neoclassical
ax = fig.add_subplot(224)

ax.set_xlabel(rtext,size=16)
ax.set_ylabel(r'Neoclassical Fluxes (MW/m$^2$)')

a=ax.plot(tgyro1.get_r(), tgyro1.get_flux_i_neo()*qgb,'b')
b=ax.plot(tgyro1.get_r(), tgyro1.get_flux_e_neo()*qgb,'r')
ax.legend((a,b),(r'$Q_D/Q_{\rm GB}$',
                 r'$Q_e/Q_{\rm GB}$'),'upper left',prop=prop)

print tgyro1.profile.keys()

plt.show()
