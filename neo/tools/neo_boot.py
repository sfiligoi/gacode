import os
import numpy as np
import sys
import string

workdir='bdir'

if len(sys.argv) < 2:
   print "python neo_boot.py <tau>"
   sys.exit()

tau = sys.argv[1]

# Prepare simulation directory
os.system('rm -rf '+workdir)
os.system('neo -g imp_c4 ; mv imp_c4 '+workdir)

# Open input.neo, append parameters, close
neoin = open(workdir+'/input.neo','a') 
neoin.write('EQUILIBRIUM_MODEL=2\n')
neoin.write('TEMP_1='+tau+'\n')
neoin.close()

# Run NEO
os.system('neo -e '+workdir)

# Harvest output
neoout = np.loadtxt(workdir+'/out.neo.transport') 
jneo=neoout[2]

neoout = open(workdir+'/out.neo.diagnostic_geo','r') 
line = neoout.readlines()[0]
ipsi=string.splitfields(line,'=')[1].rstrip()

neoout = np.loadtxt(workdir+'/out.neo.equil') 
rhostar=neoout[3]

neoout = np.loadtxt(workdir+'/out.neo.theory') 
jsauter=neoout[10]

print 'jneo     '+str(jneo)
print 'jsauter  '+str(jsauter)
print 'ipsi     '+str(float(ipsi))
print 'rhostar  '+str(rhostar)

