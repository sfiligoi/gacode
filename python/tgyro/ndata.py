import numpy as np

class tgyrodata:
    """TGYRO output data class"""

    #---------------------------------------------------------------------------#
    # Methods

    def __init__(self, sim_directory):
        """Constructor reads data from sim_directory"""
        
        self.loc_n_ion    = 0
        self.n_iterations = 0
        self.n_fields     = 0
        self.n_radial     = 0
        self.data         = {}
    
        self.dir = sim_directory
        self.getcontrol()
        #self.getresidual()
        self.getdata()

    #---------------------------------------------------------------------------#

    def getdata(self):
        """Read all tgyro-format datafiles"""
        
        self.n_ion = int(self.get_tag_value("TGYRO_N_ION"))

        self.fileparser('out.tgyro.alpha')
        self.fileparser('out.tgyro.evo_er')
        self.fileparser('out.tgyro.evo_te')
        self.fileparser('out.tgyro.evo_ti')
        self.fileparser('out.tgyro.geometry.1')
        self.fileparser('out.tgyro.geometry.2')
        self.fileparser('out.tgyro.gyrobohm')
        self.fileparser('out.tgyro.nu_rho')
        self.fileparser('out.tgyro.power_e')
        self.fileparser('out.tgyro.power_i')
        self.fileparser('out.tgyro.profile')
        self.fileparser('out.tgyro.residual')

        # Species series 
        self.fileparser('out.tgyro.evo_ne')
        self.fileparser('out.tgyro.flux_e')
        self.fileparser('out.tgyro.profile_e')
        for i in range(self.n_ion):
            for root in ['evo_n','flux_i','profile_i']:
                self.fileparser('out.tgyro.'+root+str(i+1))

    #---------------------------------------------------------------------------#

    def get_tag_value(self, tag):
        """
        Return the specified variable from input.tgyro.gen.
        tag = input.tgyro tag
        """

        datafile = file(self.dir+'/input.tgyro.gen','r')

        for line in datafile:
            try:
                if line.split()[1] == tag:
                    return float(line.split()[0])
            except IndexError:
                print "Cannot find specified input parameter: ", tag
                return 0

    #---------------------------------------------------------------------------#

    def getcontrol(self):
        """
        Read control.out to set resolutions.
        """

        import numpy as np
    
        data = np.loadtxt(self.dir+'/out.tgyro.control')

        self.n_r = int(data[0])
        self.n_evolve     = int(data[1])
        self.n_iterations = int(data[2])

    #---------------------------------------------------------------------------#

    def getresidual(self):
        """Read out.tgyro.residual
        """
        import string
        import numpy as np

        fn = 'out.tgyro.residual'
        data = open(self.dir+'/'+fn,'r').readlines()
        
        # Data dimensions 
        nr = self.n_r
        nb = self.n_iterations+1
        # 11 = 1+2*n_evolve, where n_evolve=5 (ti,te,ne,er,he)
        nc = 11
        
        numdata = np.zeros((nc,nb,nr-1),dtype=float)
        
        for ib in range(nb):
            try:
                tags=string.split(data[ib*nr]) # Contains overall residual
            except:
                print "WARNING: (data.py) out.tgyro.residual shorter than expected."
                return 0
               
            for ir in range(nr-1):
                row=string.split(data[ib*nr+ir+1])
                for ic in range(nc):
                    numdata[ic,ib,ir] = row[ic]

        self.data['residual'] = numdata
        
    #---------------------------------------------------------------------------#

    def fileparser(self,file):
        """Generic parser for standard TGYRO iteration-based file format"""

        import string
        import numpy as np

        data = open(self.dir+'/'+file,'r').readlines()

        # Data dimensions
        nr = self.n_r+2
        nc = len(string.split(data[0]))
        nb = self.n_iterations+1

        print nc

        numdata = np.zeros((nc,nb,self.n_r))
        for ib in range(nb):
            try:
                tags = string.split(data[ib*nr])
            except:
                print "WARNING: (data.py) "+file+" shorter than expected."
                return 0

            for ir in range(self.n_r):
                row = string.split(data[ib*nr+ir+2])
                numdata[:,ib,ir] = row[:]
                print row[:]

        # Populate data list
        for ic in range(nc):
            self.data[tags[ic]] = numdata[ic,:,:]
        
        print 'INFO: (data.py) Read data in '+file

