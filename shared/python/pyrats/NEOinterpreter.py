from pyrats.data import NEOData
import sys
import matplotlib.pyplot as plt

sim1 = NEOData(sys.argv[1])
sim1.scatter('SIM_2GAMMA')
sim1.scatter('SIM_2Q')
sim1.scatter('SIM_1upar')
sim1.scatter('SIM_1jboot')
plt.show()
