import os
import numpy as np
import sys

class cgyrodata:

    """CGYRO output data class."""

    def __init__(self, sim_directory):

        """Constructor reads in all simulation data."""

        self.dir = sim_directory
        self.getdata()
   

    def getdata(self):

        """Initialize data objects"""

        import numpy as np

        self.t      = np.loadtxt(self.dir+'out.cgyro.time')
        self.n_time = len(self.t)   
        print "INFO: (data.py) Read time vector in out.cgyro.time."

        #-----------------------------------------------------------------
        # Grid data is packed.  We need to unpack into sensible bits
        #
        data = np.loadtxt(self.dir+'out.cgyro.grids')

        self.n_n = int(data[0])
        self.n_species = int(data[1])
        self.n_radial  = int(data[2])
        self.n_theta   = int(data[3])
        self.n_energy  = int(data[4])
        self.n_xi      = int(data[5])
        self.m_box     = int(data[6])

        self.p = np.array(data[7:7+self.n_radial],dtype=int)
        
        mark  = 7+self.n_radial
        self.theta = np.array(data[mark:mark+self.n_theta])
        
        mark   = mark+self.n_theta
        self.energy = np.array(data[mark:mark+self.n_energy])

        mark = mark+self.n_energy
        self.xi   = np.array(data[mark:mark+self.n_xi])

        mark = mark+self.n_xi
        self.thetab = np.array(data[mark:mark+self.n_theta*self.n_radial/self.m_box])        
        print "INFO: (data.py) Read grid data in out.cgyro.grids."
        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        # Particle and energy fluxes
        #
        data = np.loadtxt(self.dir+'out.cgyro.flux_n')
        self.flux_n = np.reshape(data,(self.n_n,self.n_time),'F')
        
        data = np.loadtxt(self.dir+'out.cgyro.flux_e')
        self.flux_e = np.reshape(data,(self.n_n,self.n_time),'F')
