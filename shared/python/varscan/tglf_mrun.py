import numpy as np
import os
import string

var1 = 'RLTS_1'
vec1 = [1.0,1.1]

var2 = 'RLNS_1'
vec2 = [1.0,1.1]

outfile = open('out.tglf.array','w')
for i1 in vec1:
    for i2 in vec2:

        a = var1+'='+str(i1)
        b = var2+'='+str(i2)
        os.system('cp out.tglf.localdump input.tglf')
        os.system('echo '+a+' >> input.tglf')
        os.system('echo '+b+' >> input.tglf')
        os.system('tglf -e . > out')

        p = 0
        for line in open('out','r').readlines():
            p = p+1
            if p == 3:
                # Electron line
                elec=line.split('elec')[1].rstrip()
            if p == 4:
                # ion1 line
                ion1=line.split('ion1')[1].rstrip()

        print a+' ; '+b
        outfile.write(elec+ion1)
