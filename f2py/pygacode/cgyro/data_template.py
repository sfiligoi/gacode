import os
import numpy as np
from gacodefuncs import *
from cgyro.data import cgyrodata

sim = cgyrodata('reg01/')

# all relevant quantities available in gacode/f2py/pygacode/cgyro/data.py

print('kappa',sim.kappa)
print('ky*rho',sim.ky0)

print('omega',sim.freq[0,0,-1])
print('gamma',sim.freq[1,0,-1])
