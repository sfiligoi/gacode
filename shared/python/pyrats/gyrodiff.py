from pyrats.data import GYROData
import numpy as np
import matplotlib.pyplot as plt

sim1 = GYROData('Insert dir1 Here')
sim2 = GYROData('Insert dir2 Here')

sim1.make_diff()
sim2.make_diff()
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(sim1.diff['add dimensional specification'], 'b')
ax.plot(sim2.diff['add dimensional specification'], 'r')
ax.set_title('title')
ax.set_ylabel('ylabel')
ax.set_xlabel('xlabel')
plt.show()
