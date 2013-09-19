# Usage:
#
#   python gacode_pfile_tool.py <datafile>
#
#   NOTES: 

import sys
import numpy as np

# <datafile>
try:
    infile=sys.argv[1]
except:
    print 'Usage: python gacode_pfile_tool.py <datafile>'
    sys.exit()


outfile = open('pfile.ne','w')
for line in open(infile,'r').readlines():
    if 'psinorm' in line:
        # split() with no argument splits on whitespace
        x   = line.split()
        # Number of points
        n   = x[0]
        var = x[2]
        nvar = var.split('(')[0]
        outfile.close()
        title = 'pfile.'+nvar
        outfile = open(title,'w')
        print 'INFO: (gacode_pfile_tool.py) Extracted '+title+'.'
        outfile.write(n+'\n')
    else:
        outfile.write(line)


