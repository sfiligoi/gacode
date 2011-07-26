from pyrats.data import TGYROData
import numpy as np
import matplotlib.pyplot as plt

sim1 = TGYROData('Insert dir1 Here')
sim2 = TRYROData('Insert dir2 Here')

fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(sim1.get_something())
ax.plot(sim2.get_something())
plt.show()
