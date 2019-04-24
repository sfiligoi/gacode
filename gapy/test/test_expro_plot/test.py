import gapy
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(8,5)) ; ax = fig.add_subplot(111)
ax.set_xlabel('rho')
ax.set_ylabel('Bunit')

gapy.expro.expro_read('input.gacode')

x = gapy.expro.expro_rho
y = gapy.expro.expro_sdelta

ax.plot(x,y)

gapy.expro.expro_read('input.gacode.g')

x = gapy.expro.expro_rho
y = gapy.expro.expro_sdelta

ax.plot(x,y)

plt.show()

