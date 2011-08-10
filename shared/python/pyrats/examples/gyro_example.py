from pyrats.gyro.data import GYROData
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager
prop = matplotlib.font_manager.FontProperties(size=16)

sim1 = GYROData('/home/candy/mysim/pyrats/gyro1')
sim2 = GYROData('/home/candy/mysim/pyrats/gyro2')

fig = plt.figure(figsize=(12,5))
fig.subplots_adjust(left=0.16,right=0.94)
fig.subplots_adjust(bottom=0.12,top=0.92)

#-----------------------------------------------------------
ax = fig.add_subplot(121)

#gbflux: n_kinetic x n_field x 4 x n_time
 
# ions
a=ax.plot(sim1.t['(cbar_s/a)t'],sim1.gbflux[0, 0, 1, :], 'k')
b=ax.plot(sim2.t['(cbar_s/a)t'],sim2.gbflux[0, 0, 1, :], 'k--')
# electrons
c=ax.plot(sim1.t['(cbar_s/a)t'],sim1.gbflux[1, 0, 1, :], 'b')
d=ax.plot(sim2.t['(cbar_s/a)t'],sim2.gbflux[1, 0, 1, :], 'b--')
ax.set_xlabel(r'$(c_s/a)t$',size=16)
ax.legend((a,b,c,d),(r'$Q_i^{(1)}$',
                     r'$Q_i^{(2)}$',
                     r'$Q_e^{(1)}$',
                     r'$Q_e^{(2)}$'),prop=prop)
#-----------------------------------------------------------

#-----------------------------------------------------------
ax = fig.add_subplot(122)

#gbflux_n: n_kinetic x n_field x 4 x n_n x n_time

n_n = int(sim1.profile['n_n'])

# Sim-1
n_time = int(sim1.t['n_time'])
t1 = n_time/10
t2 = n_time

f = np.zeros(n_n)
# ions
for i in range(n_n):
    f[i] = sim1.gbflux_n[0,0,1,i,t1:t2].sum()/(t2-t1)
a=ax.plot(sim1.profile['kt_rho'],f,'k')

# electrons
for i in range(n_n):
    f[i] = sim1.gbflux_n[1,0,1,i,t1:t2].sum()/(t2-t1)
c=ax.plot(sim1.profile['kt_rho'],f,'b')

# Sim-2
n_time = int(sim2.t['n_time'])
t1 = n_time/10
t2 = n_time

f = np.zeros(n_n)
# ions
for i in range(n_n):
    f[i] = sim2.gbflux_n[0,0,1,i,t1:t2].sum()/(t2-t1)
b=ax.plot(sim2.profile['kt_rho'],f,'k')

# electrons
for i in range(n_n):
    f[i] = sim2.gbflux_n[1,0,1,i,t1:t2].sum()/(t2-t1)
d=ax.plot(sim2.profile['kt_rho'],f,'b')

ax.set_xlabel(r'$k_\theta \rho_s$',size=16)
ax.legend((a,b,c,d),(r'$Q_i^{(1)}$',
                     r'$Q_i^{(2)}$',
                     r'$Q_e^{(1)}$',
                     r'$Q_e^{(2)}$'),prop=prop)
#-----------------------------------------------------------

plt.show()
