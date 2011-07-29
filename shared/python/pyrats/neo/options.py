from pyrats.neo.data import NEOData
import sys
import numpy as np

sim1 = NEOData(sys.argv[1])
keys = []

for k1, v1 in sim1.transport.iteritems():
    s = 0
    for k2, v2 in v1.data.iteritems():
        for item in v2:
            s = s + np.sum(item)
    if s != 0:
#Unless every entry in v2 is zero, it gets added to the list.
        keys.append(k1)
keys.sort()
for k in keys:
    print k
