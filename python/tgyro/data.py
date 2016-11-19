import os
import sys
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
        self.n_ion = int(self.get_tag_value("LOC_N_ION"))
        self.getcontrol()

        # Manage old data
        if os.path.isfile(self.dir+'/out.tgyro.profile_i1'):
            print "INFO: (data.py) Detected new multi-ion version of TGYRO [GOOD]."
            self.getdata()
        else:
            print "INFO: (data.py) Detected old static ion version of TGYRO."
            print "      --------> Using getolddata function ..."
            self.getolddata() 
        
    #---------------------------------------------------------------------------#

    def getolddata(self):
        """Reconstruct new datafiles from old dataset"""

        self.fileparser('out.tgyro.alpha')
        self.fileparser('out.tgyro.geometry.1')
        self.fileparser('out.tgyro.geometry.2')
        self.fileparser('out.tgyro.gyrobohm')
        self.fileparser('out.tgyro.nu_rho',onerow=True)
        self.fileparser('out.tgyro.power_e')
        self.fileparser('out.tgyro.power_i')

        self.fileparser('out.tgyro.flux_e')
        self.fileparser('out.tgyro.flux_i')
        self.fileparser('out.tgyro.flux_target')
        self.fileparser('out.tgyro.gradient')
        self.fileparser('out.tgyro.mflux_e')
        self.fileparser('out.tgyro.mflux_i')
        self.fileparser('out.tgyro.profile')
        self.data['pflux_i1_tot'] = self.data['pflux_i_neo']+self.data['pflux_i_tur']
        self.data['pflux_i1_target'] = 0.0*self.data['pflux_i_neo']
        for i in np.arange(1,self.n_ion):
            self.data['pflux_i'+str(i)+'_tot'] = 0.0
        self.data['a/Lte']  = self.data['a/LTe']
        self.data['a/Lti1'] = self.data['a/LTi']
        self.data['ti1'] = self.data['ti']
        self.data['a/Lni1'] = self.data['a/Lni']
        self.data['ni1'] = self.data['ni']

    def getdata(self):
        """Read all tgyro-format datafiles"""
        
        self.fileparser('out.tgyro.alpha')
        self.fileparser('out.tgyro.geometry.1')
        self.fileparser('out.tgyro.geometry.2')
        self.fileparser('out.tgyro.gyrobohm')
        self.fileparser('out.tgyro.nu_rho')
        self.fileparser('out.tgyro.power_e')
        self.fileparser('out.tgyro.power_i')

        self.fileparser('out.tgyro.residual')
        self.fileparser('out.tgyro.evo_er')
        self.fileparser('out.tgyro.evo_te')
        self.fileparser('out.tgyro.evo_ti')
        self.fileparser('out.tgyro.profile')

        # Species series 
        self.fileparser('out.tgyro.evo_ne')
        self.fileparser('out.tgyro.flux_e')
        self.fileparser('out.tgyro.profile_e')
        for i in range(self.n_ion):
            for root in ['evo_n','flux_i','profile_i']:
                self.fileparser('out.tgyro.'+root+str(i+1))

        print self.data.keys()

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

    def fileparser(self,file,onerow=False):
        """Generic parser for standard TGYRO iteration-based file format"""

        import string
        import numpy as np

        data = open(self.dir+'/'+file,'r').readlines()

        # Data dimensions
        if onerow == False:
            ix = 2
        else:
            ix = 1
        nr = self.n_r+ix
        nc = len(string.split(data[0]))
        nb = self.n_iterations+1

        numdata = np.zeros((nc,nb,self.n_r))
        for ib in range(nb):
            try:
                tags = string.split(data[ib*nr])
            except:
                print "WARNING: (data.py) "+file+" shorter than expected."
                return 0

            for ir in range(self.n_r):
                row = string.split(data[ib*nr+ir+ix])
                numdata[:,ib,ir] = row[:]

        # Populate data list
        for ic in range(nc):
            self.data[tags[ic]] = numdata[ic,:,:]
        
        print 'INFO: (data.py) Read data in '+file

