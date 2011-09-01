class NEOData:
    """A class for management of NEO output data.

    Data:

    Example Usage:
        >>> from pyrats.neo.data import NEOData
        >>> sim1 = NEOData('example_directory')
    """
    
    #---------------------------------------------------------------------------#
    # Methods

    def __init__(self, sim_directory):
        """Constructor reads in data from sim_directory and creates new object.
        """

        self.init_data()
        self.set_directory(sim_directory)
        self.read_grid()
        self.read_equil()
        self.read_theory()
        self.read_transport()
        self.read_transport_exp()
        self.read_vel()

    #---------------------------------------------------------------------------#

    def init_data(self):
        """Initialize object data."""

        self.dirname       = ""
        self.grid          = {}
        self.equil         = {}
        self.theory        = {}
        self.transport     = {}
        self.transport_exp = {}
        self.vel           = []

    #---------------------------------------------------------------------------#

    def set_directory(self, path):
        """Set the simulation directory."""

        from os.path import expanduser, expandvars
        self.dirname = expanduser(expandvars(path))
    
    #---------------------------------------------------------------------------#

    def read_grid(self):
        """Reads out.neo.grid"""

        import sys
        import numpy as np
 
        try:
            data = np.loadtxt(self.dirname+'/out.neo.grid')
        except:
            print "ERROR (NEOData): Fatal error!  Missing out.neo.grid."
            sys.exit()
                
        self.grid['n_species'] = int(data[0])
        self.grid['n_energy']  = int(data[1])
        self.grid['n_xi']      = int(data[2])
        self.grid['n_theta']   = int(data[3])
        self.grid['theta']     = data[4:4+int(data[3])]
        self.grid['n_radial']  = int(data[4+int(data[3])])
        self.grid['r_over_a']  = data[5+int(data[3]):]

    #---------------------------------------------------------------------------#

    def read_equil(self):
        """Read out.neo.equil."""

        import sys
        import numpy as np
 
        try:
            equil = np.loadtxt(self.dirname+'/out.neo.equil')
        except:
            print "ERROR (NEOData): Fatal error!  Missing out.neo.equil."
            sys.exit()
            
        n_spec = self.grid['n_species']

        self.equil['r_over_a']      = equil[:,0]
        self.equil['dphidr']        = equil[:,1]
        self.equil['q']             = equil[:,2]
        self.equil['rho_star']      = equil[:,3]
        self.equil['R0_over_a']     = equil[:,4]
        self.equil['omega0']        = equil[:,5]
        self.equil['domega0dr']     = equil[:,6]
        self.equil['n_norm']        = equil[:,7]
        self.equil['T_norm']        = equil[:,8]
        self.equil['v_norm_over_a'] = equil[:,9]
        self.equil['n']             = equil[:,10+0*n_spec:10+1*n_spec]
        self.equil['T']             = equil[:,10+1*n_spec:10+2*n_spec]
        self.equil['a_over_ln']     = equil[:,10+2*n_spec:10+3*n_spec]
        self.equil['a_over_lt']     = equil[:,10+3*n_spec:10+4*n_spec]
        self.equil['tauinv_norm']   = equil[:,10+4*n_spec:10+5*n_spec]

    #---------------------------------------------------------------------------#

    def read_theory(self):
        """Reads out.neo.theory."""

        import sys
        import numpy as np
 
        try:
            data = np.loadtxt(self.dirname+'/out.neo.theory')
        except:
            print "ERROR (NEOData): Fatal error!  Missing out.neo.theory."
            sys.exit()
            
        n_spec = self.grid['n_species']

        self.theory['r_over_a'] = data[:,0]
        self.theory['HHGamma']  = data[:,1]
        self.theory['HHQi']     = data[:,2]
        self.theory['HHQe']     = data[:,3]
        self.theory['HHjparB']  = data[:,4]
        self.theory['HHk']      = data[:,5]
        self.theory['HHuparB']  = data[:,6]
        self.theory['HHvtheta'] = data[:,7]
        self.theory['CHQi']     = data[:,8]
        self.theory['TGQi']     = data[:,9]
        self.theory['SjparB']   = data[:,10]
        self.theory['Sk']       = data[:,11]
        self.theory['SuparB']   = data[:,12]
        self.theory['Svtheta']  = data[:,13]
        self.theory['HRphisq']  = data[:,14]
        self.theory['HSGamma']  = data[:,15+0*n_spec:15+1*n_spec]
        self.theory['HSQ']      = data[:,15+1*n_spec:15+2*n_spec]

    #---------------------------------------------------------------------------#

    def read_transport(self):
        """Reads out.neo.transport."""
        
        import sys
        import numpy as np
 
        try:
            data = np.loadtxt(self.dirname+'/out.neo.transport')
        except:
            print "ERROR (NEOData): Fatal error!  Missing out.neo.transport."
            sys.exit()
            
        n_spec = self.grid['n_species']

        self.transport['r_over_a'] = data[:,0]
        self.transport['phisq']    = data[:,1]
        self.transport['jparB']    = data[:,2]
        self.transport['vtheta']   = data[:,3]
        self.transport['uparB']    = data[:,4]
        self.transport['Gamma']    = data[:,5+0*n_spec:5+1*n_spec]
        self.transport['Q']        = data[:,5+1*n_spec:5+2*n_spec]
        self.transport['Pi']       = data[:,5+2*n_spec:5+3*n_spec]
        self.transport['uparB']    = data[:,5+3*n_spec:5+4*n_spec]
        self.transport['k']        = data[:,5+4*n_spec:5+5*n_spec]
        self.transport['K']        = data[:,5+5*n_spec:5+6*n_spec]
        self.transport['vtheta']   = data[:,5+6*n_spec:5+7*n_spec]
        self.transport['vphi']     = data[:,5+7*n_spec:5+8*n_spec]

    #---------------------------------------------------------------------------#

    def read_transport_exp(self):
        """Reads out.neo.transport.exp."""
        
        import sys
        import numpy as np
 
        try:
            data = np.loadtxt(self.dirname+'/out.neo.transport_exp')
        except:
            print "ERROR (NEOData): Fatal error!  Missing out.neo.transport.exp."
            sys.exit()
            
        n_spec = self.grid['n_species']

        self.transport_exp['r']      = data[:,0]
        self.transport_exp['phisq']  = data[:,1]
        self.transport_exp['jparB']  = data[:,2]
        self.transport_exp['vtheta'] = data[:,3]
        self.transport_exp['uparB']  = data[:,4]
        self.transport_exp['Gamma']  = data[:,5+0*n_spec:5+1*n_spec]
        self.transport_exp['Q']      = data[:,5+1*n_spec:5+2*n_spec]
        self.transport_exp['Pi']     = data[:,5+2*n_spec:5+3*n_spec]
        self.transport_exp['uparB']  = data[:,5+3*n_spec:5+4*n_spec]
        self.transport_exp['k']      = data[:,5+4*n_spec:5+5*n_spec]
        self.transport_exp['K']      = data[:,5+5*n_spec:5+6*n_spec]
        self.transport_exp['vtheta'] = data[:,5+6*n_spec:5+7*n_spec]
        self.transport_exp['vphi']   = data[:,5+7*n_spec:5+8*n_spec]

    #---------------------------------------------------------------------------#

    def read_vel(self):
        """Reads out.neo.vel."""
        
        import sys
        import numpy as np
 
        try:
            data = np.loadtxt(self.dirname+'/out.neo.vel')
        except:
            print "ERROR (NEOData): Missing out.neo.vel."
            sys.exit()
                   
        self.vel = data.reshape((self.grid['n_radial'],
                                 self.grid['n_species'],
                                 self.grid['n_theta']))
