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

        """Initialize smaller data objects (don't load larger ones)"""

        import numpy as np
        import time

        #-----------------------------------------------------------------
        # Read time vector.
        #
        data = np.loadtxt(self.dir+'out.cgyro.time')
        self.t      = data[:,0]
        self.error  = data[:,1]
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
        self.n_theta_plot = int(data[7])
        self.m_box     = int(data[8])
        self.length    = float(data[9])
        # Set l to last data index plus one.
        l=10

        self.p = np.array(data[l:l+self.n_radial],dtype=int)
        self.kx = 2*np.pi*self.p/self.length
        
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
            data = np.fromfile(self.dir+'out.cgyro.freq',dtype='float',sep=" ")
            self.freq = np.reshape(data,(2,self.n_n,nt),'F')
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
                                       self.n_species,self.n_xi,self.n_energy,nt),'F')
            self.hb = self.hb/np.max(self.hb)
            print "INFO: (data.py) Read data in out.cgyro.hb."
        except:
            pass
        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        # Particle and energy fluxes
        #
        try:
            start = time.time()
            data = np.fromfile(self.dir+'out.cgyro.kxky_flux_n',dtype='float',sep=" ")
            end = time.time()
            self.flux_n = np.reshape(data,(self.n_radial,self.n_species,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.kxky_flux_n. TIME = "+str(end-start)
        except:
            pass

        try:
            start = time.time()
            data = np.fromfile(self.dir+'out.cgyro.kxky_flux_e',dtype='float',sep=" ")
            end = time.time()
            self.flux_e = np.reshape(data,(self.n_radial,self.n_species,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.kxky_flux_e. TIME = "+str(end-start)
        except:
            pass
        #-----------------------------------------------------------------


    def getbig(self):

        """Larger files"""

        import numpy as np
        import time
        
        # Convenience definition
        nt = self.n_time

        #-----------------------------------------------------------------
        # Read complex fields
        #
        try:
            start = time.time()
            data = np.fromfile(self.dir+'out.cgyro.kxky_phi',dtype='float',sep=" ")
            end = time.time()
            self.phi = np.reshape(data,(2,self.n_radial,self.n_n,nt),'F')
            self.phisq = self.phi[0,:,:,:]**2+self.phi[1,:,:,:]**2
            print "INFO: (data.py) Read data in out.cgyro.kxky_phi. TIME = "+str(end-start)
        except:
            pass
        
        try:
            start = time.time()
            data = np.fromfile(self.dir+'out.cgyro.kxky_n',dtype='float',sep=" ")
            end = time.time()
            self.n = np.reshape(data,(2,self.n_radial,self.n_species,self.n_n,nt),'F')
            self.nsq = self.n[0,:,:,:,:]**2+self.n[1,:,:,:,:]**2
            print "INFO: (data.py) Read data in out.cgyro.kxky_n. TIME = "+str(end-start)
        except:
            pass
        #-----------------------------------------------------------------
