import os
import sys
import numpy as np
from gacodefuncs import *

def analyze_resolution(N,B,twoM):
    """
    Perform combinatoric mixing + wrapping analysis for n=1..N.
    
    Parameters
    ----------
    N     : max n index
    B     : box size
    twoM  : radial domain size (p runs from -M..M-1)
    
    """
    
    m = twoM//2
    
    print("  n     l    score   rating ")
    print("----  ----  ------- --------")
   
    for n in range(1,N+1):
        l = B*n
        
        # -------------------------
        # Mixing score: l relative to 2M
        # -------------------------
        score = l/twoM

        if score < 0.05:
           rating = "*"*(10+int((0.05/score)))
        elif score < 0.5:
           rating = '*'*int(0.5/score)
        elif score < 1.0:
           rating = "poorly resolved"
        else:
           rating = "unresolved"
                
        
        print("{:3d}  {:4d}    {:.3f}   {}".format(n,l,score,rating))

n_radial   = int(parameter_val('input.cgyro.gen','N_RADIAL'))
n_toroidal = int(parameter_val('input.cgyro.gen','N_TOROIDAL'))
box_size   = int(parameter_val('input.cgyro.gen','BOX_SIZE'))

if n_toroidal == 1:
   n = 1
else:
   n = n_toroidal-1
   
analyze_resolution(n,box_size,n_radial)
