"""This script gives help for the requested TGYROData data object."""

import sys
from pyrats.tgyro.data import TGYROData

sim1 = TGYROData(sys.argv[1])

if sys.argv[2] == 'data':
    print sorted(sim1.data.keys())
elif sys.argv[2] == 'docstrings':
    help(TGYROData)
