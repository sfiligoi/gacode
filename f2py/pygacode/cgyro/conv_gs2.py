import sys
import os
import numpy as np
from gacodefuncs import *
from .data import cgyrodata

# read basic simulation data
sim = cgyrodata('./')

print('B_gs2/B_unit = {:f}'.format(sim.b_gs2))

c = sim.b_gs2

print('Beta_e_unit = {:f}'.format(sim.betae_unit))
print('Beta_e_gs2 = {:f}'.format(sim.betae_unit/c**2))

