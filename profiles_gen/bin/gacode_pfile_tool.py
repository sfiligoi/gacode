# Usage:
#
#   python gacode_pfile_tool.py <datafile>

import sys
import numpy as np

try:
    infile=sys.argv[1]
except:
    print('Usage: python gacode_pfile_tool.py <datafile>')
    sys.exit()

outfile = open('pfile.ne','w')
with open(infile,'r') as fin:
    in_lines = fin.readlines()
for line in in_lines:
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
        print('INFO: (gacode_pfile_tool.py) Extracted '+title+'.')
        outfile.write(n+'\n')
    elif 'SPECIES' in line:
        outfile.close()
        title = 'pfile.species'
        outfile = open(title,'w')
        print('INFO: (gacode_pfile_tool.py) Found ion species header.')
        print('INFO: (gacode_pfile_tool.py) Extracted '+title+'.')
        n = line[0]
        outfile.write(n+'\n')
    else:
        outfile.write(line)
outfile.close()
