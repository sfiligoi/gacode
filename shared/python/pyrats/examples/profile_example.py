from pyrats.profiles_gen.data import profiles_genData
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager
prop = matplotlib.font_manager.FontProperties(size=16)

path = '/home/candy/mysim/pyrats/'
pro1 = profiles_genData(path+'neo1')

fig = plt.figure(figsize=(12,8))
fig.subplots_adjust(left=0.16,right=0.94)
fig.subplots_adjust(bottom=0.12,top=0.92)

#pro1.plot('ne', fignum=1, plotcounter=1)

#-----------------------------------------------------------------------
ax = fig.add_subplot(221)

ax.set_xlabel(r'$r\,(m)$',size=16)
ax.set_ylabel(r'$10^{19}/{\rm m}^3$',size=16)

a = ax.plot(pro1.data['rmin (m)'],pro1.data['ne (10^19/m^3)'],'b')
b = ax.plot(pro1.data['rmin (m)'],pro1.data['ni_1 (10^19/m^3)'],'r')
c = ax.plot(pro1.data['rmin (m)'],np.double(pro1.data['ni_2 (10^19/m^3)'])*6.0,'g')
ax.legend((a,b,c),(r'$n_e$',r'$n_{D}$',r'$n_C Z_C$'),
 'center left',prop=prop)
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
ax = fig.add_subplot(222)

ax.set_xlabel(r'$r\,(m)$',size=16)
ax.set_ylabel(r'${\rm keV}$',size=16)

a = ax.plot(pro1.data['rmin (m)'],pro1.data['Te (keV)'],'b')
b = ax.plot(pro1.data['rmin (m)'],pro1.data['Ti_1 (keV)'],'r')
ax.legend((a,b),(r'$T_e$',r'$T_{D}=T_{C}$'),
 prop=prop)
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
ax = fig.add_subplot(223)

ax.set_xlabel(r'$r\,(m)$',size=16)

print pro1.data.keys()

rmaj = np.double(pro1.data['rmaj (m)'])

a = ax.plot(pro1.data['rmin (m)'],pro1.data['kappa (-)'],'b')
b = ax.plot(pro1.data['rmin (m)'],pro1.data['delta (-)'],'r')
c = ax.plot(pro1.data['rmin (m)'],pro1.data['zeta (-)'],'g')
d = ax.plot(pro1.data['rmin (m)'],rmaj/max(rmaj),'c')
ax.legend((a,b,c,d),(r'$\kappa$',r'$\delta$',r'$\zeta$',r'$R/a$'),
 'lower left',prop=prop)
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
ax = fig.add_subplot(224)

ax.set_xlabel(r'$r\,(m)$',size=16)
ax.set_ylabel(r'${\rm km/s}$',size=16)

v1 = np.double(pro1.data['vtor_1 (m/s)'])/1000.0
v2 = np.double(pro1.data['vtor_2 (m/s)'])/1000.0

a = ax.plot(pro1.data['rmin (m)'],v1,'b')
b = ax.plot(pro1.data['rmin (m)'],v2,'r')
ax.legend((a,b),(r'$v_{\rm \varphi,D}$',r'$v_{\rm \varphi,C}$'),
 'lower left',prop=prop)
#-----------------------------------------------------------------------

#plt.show()
plt.savefig('profile.pdf')
