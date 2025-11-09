import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import AutoMinorLocator
import matplotlib.gridspec as gridspec
from cgyro.data import cgyrodata
from sDMD.sDMD import sDMD 


DTYPE='float64'


mydir='/users/path/to/CGYRO_case/'
time = np.loadtxt(".../time.txt")
theta = np.loadtxt(".../theta.txt")
tmp_re = np.loadtxt(".../phi_re.txt")
tmp_im = np.loadtxt(".../phi_im.txt")


Nt = len(time)
Ntheta = len(theta)
tmp = np.zeros([Nt, Ntheta], dtype=DTYPE)
phi_r = tmp_re.reshape(tmp.shape)
phi_i = tmp_im.reshape(tmp.shape)

# how many points to keep in gyrokinetic observables 
k = 5 #25 #50



Y = phi_r[:,:]
Y = Y.T
y_start, y_rest = Y[:,:k], Y[:,k:]
sdmd = sDMD(y_start, rmin=2, rmax=3, f=5, s=1)

Y_predict = []
for y in y_rest.T:
    out = sdmd.update(y)

    prediction = sdmd.Uy @ sdmd.A @ sdmd.Ux.T @ y
    Y_predict.append(prediction)

Y_predict = np.array(Y_predict).T
y_actual = y_rest[:, sdmd.f:]
y_predict = Y_predict

n = len(y_predict[15,k:])



fig = plt.figure()
gs = gridspec.GridSpec(1,1)
ax1 = fig.add_subplot(gs[0,0])

ax1.plot(time, phi_r[:,15], '-', color='red', label='GKs')
ax1.plot(time[-n:], y_predict[15,k:], color='tab:blue', label='DMD')

ax1.tick_params(labelsize=15)
ax1.legend(loc=3, fontsize=12.5)
fig.tight_layout()
plt.show()



