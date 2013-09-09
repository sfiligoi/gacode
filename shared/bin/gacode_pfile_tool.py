# Usage:
#
#   python gacode_pfile_tool.py <datafile>
#
#   NOTES: 

import sys
import numpy as np
from itertools import islice

# <datafile>
try:
    infile=sys.argv[1]
except:
    print 'Usage: python gacode_pfile_tool.py <datafile>'
    sys.exit()


outfile = open('pfile.ne','w')
for line in open(infile,'r').readlines():
    if 'psinorm' in line:
        x   = line.split(' ')
        n   = x[0]
        var = x[2]
        nvar = var.split('(')[0]
        outfile.close()
        title = 'pfile.'+nvar
        outfile = open(title,'w')
        outfile.write(n+'\n')
    else:
        outfile.write(line)


