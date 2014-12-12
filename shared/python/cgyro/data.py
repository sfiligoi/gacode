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

        #-----------------------------------------------------------------
        # Read time vector.
        #
        self.t      = np.loadtxt(self.dir+'out.cgyro.time')
        self.n_time = len(self.t)   
        print
        print "INFO: (data.py) Read time vector in out.cgyro.time."
        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        # Read grid data.
        # 
        # NOTE: Grid data is packed, so unpack into sensible bits
        #
        data = np.loadtxt(self.dir+'out.cgyro.grids')

        self.n_n       = int(data[0])
        self.n_species = int(data[1])
        self.n_field   = int(data[2])
        self.n_radial  = int(data[3])
        self.n_theta   = int(data[4])
        self.n_energy  = int(data[5])
        self.n_xi      = int(data[6])
        self.m_box     = int(data[7])
        # Set l to last data index plus one.
        l=8

        self.p = np.array(data[l:l+self.n_radial],dtype=int)
        
        mark = l+self.n_radial
        self.theta = np.array(data[mark:mark+self.n_theta])
        
        mark = mark+self.n_theta
        self.energy = np.array(data[mark:mark+self.n_energy])

        mark = mark+self.n_energy
        self.xi   = np.array(data[mark:mark+self.n_xi])

        mark = mark+self.n_xi
        self.thetab = np.array(data[mark:mark+self.n_theta*self.n_radial/self.m_box])        
         
        mark = mark+self.n_theta*self.n_radial/self.m_box
        self.ky = np.array(data[mark:mark+self.n_n])

        print "INFO: (data.py) Read grid data in out.cgyro.grids."
        #-----------------------------------------------------------------

        # Convenience definition
        nt = self.n_time

        #-----------------------------------------------------------------
        # Linear frequency
        #
        try:
            data = np.loadtxt(self.dir+'out.cgyro.freq')
            self.freq = data
            print "INFO: (data.py) Read data in out.cgyro.freq."
        except:
            pass

        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        # Read ballooning potentials
        #
        try:
            data = np.loadtxt('out.cgyro.phib')
            self.phib = np.reshape(data,(2,self.n_theta*self.n_radial/self.m_box,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.phib."
        except:
            pass
        try:
            data = np.loadtxt('out.cgyro.aparb')
            self.aparb = np.reshape(data,(2,self.n_theta*self.n_radial/self.m_box,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.aparb."
        except:
            pass
        try:
            data = np.loadtxt('out.cgyro.bparb')
            self.bparb = np.reshape(data,(2,self.n_theta*self.n_radial/self.m_box,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.bparb."
        except:
            pass
        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        # Ballooning distribution
        #
        try:
            data = np.loadtxt(self.dir+'out.cgyro.hb')
            self.hb = np.reshape(data,(2,self.n_radial*self.n_theta/self.m_box,
                                       self.n_species,self.n_xi,self.n_energy),'F')
            self.hb = self.hb/np.max(self.hb)
            print "INFO: (data.py) Read data in out.cgyro.hb."
        except:
            pass
        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        # Read standard potentials
        #
        try:
            data = np.loadtxt(self.dir+'out.cgyro.phi')
            self.phi = np.reshape(data,(2,self.n_theta,self.n_radial,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.phi."
        except:
            pass
        try:
            data = np.loadtxt(self.dir+'out.cgyro.apar')
            self.apar = np.reshape(data,(2,self.n_theta,self.n_radial,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.apar."
        except:
            pass
        try:
            data = np.loadtxt(self.dir+'out.cgyro.bpar')
            self.bpar = np.reshape(data,(2,self.n_theta,self.n_radial,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.bpar."
        except:
            pass
        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        # Read powers
        #
        try:
            data = np.loadtxt(self.dir+'out.cgyro.pwr_phi')
            self.pwr_phi = np.reshape(data,(self.n_radial,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.pwr_phi."
        except:
            pass
        try:
            data = np.loadtxt(self.dir+'out.cgyro.pwr_apar')
            self.pwr_apar = np.reshape(data,(self.n_radial,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.pwr_apar."
        except:
            pass
        try:
            data = np.loadtxt(self.dir+'out.cgyro.pwr_bpar')
            self.pwr_bpar = np.reshape(data,(self.n_radial,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.pwr_bpar."
        except:
            pass
        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        # Particle and energy fluxes
        #
        try:
            data = np.loadtxt(self.dir+'out.cgyro.flux_n')
            self.flux_n = np.reshape(data,(self.n_species,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.flux_n."
        except:
            pass

        try:
            data = np.loadtxt(self.dir+'out.cgyro.flux_e')
            self.flux_e = np.reshape(data,(self.n_species,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.flux_e."
        except:
            pass
      #-----------------------------------------------------------------


