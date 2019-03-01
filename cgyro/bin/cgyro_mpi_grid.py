import os
import sys
import numpy as np

if len(sys.argv) == 1:
   print('USAGE: python cgyro_mpi_grid.py')
   print('       (assumes input.cgyro.gen in CWD')
   sys.exit()

def gcd(m,n):

   a = m ; b = n

   if a < b:
      c = a
      a = b
      b = c

  do          
     c = mod(a, b)    
     if (c == 0) exit
     a = b         
     b = c 
  enddo

  d = b

  return d

n_xi=16
n_energy=8
n_species=2
n_radial=64
n_theta=24

# Velocity-space (v) and configuration-space (c) dimensions
nv = n_energy*n_xi*n_species
nc = n_radial*n_theta

d = gcd(nv,nc)

#-------------------------------------------------------------------------
# MPI diagnostics need to come early
#
print('Parallelization and distribution diagnostics')
#write(io,'(a,i5)') '         nv: ',nv
#write(io,'(a,i5)') '         nc: ',nc
#write(io,*)
#write(io,*) '          [coll]     [str]      [NL] '
#write(io,*) ' n_MPI    nc_loc    nv_loc   n_split '
#write(io,*) '------    ------    ------   ------- '
for it in range(d*n_toroidal):
   j = it+1
   if np.mod(d*n_toroidal,j) == 0 .and. np.mod(j,n_toroidal) == 0:
      n_proc_1 = j//n_toroidal
      nc_loc = nc//n_proc_1           
      nv_loc = nv//n_proc_1           
      nsplit = 1+(nv_loc*n_theta-1)/n_toroidal
      print(j,nc_loc,nv_loc,nsplit)


