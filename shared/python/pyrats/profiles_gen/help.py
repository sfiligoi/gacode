"""This script gives help for the requested profiles_genData data object."""

import sys
from pyrats.profiles_gen.data import profiles_genData

sim1 = profiles_genData(sys.argv[1])

if sys.argv[2] == 'data':
    print sorted(sim1.data.keys())
elif sys.argv[2] == 'docstrings':
    help(profiles_genData)
