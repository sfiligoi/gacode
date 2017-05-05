import os
import numpy as np
import sys
import string

workdir='bdir'

if len(sys.argv) < 3:
   print "python neo_boot.py <tau> <rmin>"
   sys.exit()

tau  = sys.argv[1]
rmin = sys.argv[2]

# Prepare simulation directory
os.system('rm -rf '+workdir)
os.system('neo -g imp_c4 ; mv imp_c4 '+workdir)

list = ['DLNNDR_1',
        'DLNNDR_2',
        'DLNNDR_3',
        'DLNTDR_1',
        'DLNTDR_2',
        'DLNTDR_3']

# Open input.neo, append parameters, close
neoin = open(workdir+'/input.neo','a') 
neoin.write('EQUILIBRIUM_MODEL=2\n')

# Set input: tau
neoin.write('TEMP_1='+tau+'\n')
neoin.write('TEMP_3='+tau+'\n')

# Set input: r_min
neoin.write('RMIN_OVER_A='+rmin+'\n')
neoin.close()

c = []

# Run NEO
for i in range(6):
   neoin = open(workdir+'/input.neo','a') 
   # Overlay gradients
   neoin.write('# '+str(i)+'\n')
   for j in range(6):
      if j == i:
         neoin.write(list[j]+'=1\n')
      else:        
         neoin.write(list[j]+'=0\n')
   neoin.close()
   os.system('neo -e '+workdir)

   # Harvest output
   neoout = np.loadtxt(workdir+'/out.neo.transport') 
   jneo=neoout[2]

   neoout = open(workdir+'/out.neo.diagnostic_geo','r') 
   line = neoout.readlines()[0]
   ipsi = float(string.splitfields(line,'=')[1].rstrip())

   neoout = np.loadtxt(workdir+'/out.neo.equil') 
   rhostar=neoout[3]

   neoout = np.loadtxt(workdir+'/out.neo.theory') 
   jsauter=neoout[10]

   print 'jneo        '+str(jneo)
   print 'jsauter     '+str(jsauter)
   print 'I*Psi*rho_* '+str(ipsi*rhostar)

   c.append(jneo/(ipsi*rhostar))

print list
print c
