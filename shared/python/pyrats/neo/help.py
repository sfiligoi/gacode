"""This script gives help for the requested NEOData data object."""

import sys
from pyrats.neo.data import NEOData

sim1 = NEOData(sys.argv[1])

if sys.argv[2] == 'transport':
    print sorted(sim1.transport.keys())
elif sys.argv[2] == 'HH_theory':
    print sorted(sim1.HH_theory.keys())
elif sys.argv[2] == 'CH_theory':
    print sorted(sim1.CH_theory.keys())
elif sys.argv[2] == 'TG_theory':
    print sorted(sim1.TG_theory.keys())
elif sys.argv[2] == 'S_theory':
    print sorted(sim1.S_theory.keys())
elif sys.argv[2] == 'HR_theory':
    print sorted(sim1.HR_theory.keys())
elif sys.argv[2] == 'HS_theory':
    print sorted(sim1.HS_theory.keys())
elif sys.argv[2] == 'docstrings':
    help(NEOData)
