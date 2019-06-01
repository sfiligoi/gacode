import gapy
import numpy as np
import matplotlib.pyplot as plt

gapy.expro.expro_read('input.gacode')
nexp = int(gapy.expro.expro_n_exp) 

fig = plt.figure(figsize=(8,5)) ; ax = fig.add_subplot(111)
ax.set_xlabel('rho')
ax.set_ylabel('Bunit')

x = gapy.expro.expro_rho
y = gapy.expro.expro_bunit

ax.plot(x,y)

plt.show()

