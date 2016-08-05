import numpy as np

class tgyrodata:
    """TGYRO output data class."""

    #---------------------------------------------------------------------------#
    # Methods

    def __init__(self, sim_directory):
        """Constructor reads in data from sim_directory and creates new object.
        """
        
        self.loc_n_ion    = 0
        self.n_iterations = 0
        self.n_fields     = 0
        self.n_radial     = 0
        self.data         = {}
    
        self.dir = sim_directory
        self.getdata()

    #---------------------------------------------------------------------------#

    def getdata(self):
        """Read in object data."""
        
        self.loc_n_ion = int(self.get_tag_value("LOC_N_ION"))
        self.read_control()
        self.fileparser('out.tgyro.gyrobohm')
        self.fileparser('out.tgyro.flux_e')
        self.fileparser('out.tgyro.flux_i')
        self.fileparser('out.tgyro.flux_target')
        self.fileparser('out.tgyro.mflux_e')
        self.fileparser('out.tgyro.mflux_i')
        self.fileparser('out.tgyro.profile')
        self.fileparser('out.tgyro.gradient')
        self.fileparser('out.tgyro.geometry.1')
        self.fileparser('out.tgyro.geometry.2')
        self.fileparser('out.tgyro.power_e')
        self.fileparser('out.tgyro.power_i')
        for i in range(2,self.loc_n_ion+1):
            si = '%d'%i
            for fn in ['profile','mflux_i','flux_i']:
                self.fileparser('out.tgyro.'+fn+si,spec_num=i)
        self.get_residual()
    #---------------------------------------------------------------------------#

    def get_tag_value(self, tag):
        """
        Return the specified variable from input.tgyro.gen.
        tag = input.tgyro tag
        """

        datafile = file(self.dirname+'/input.tgyro.gen','r')

        for line in datafile:
            try:
                if line.split()[1] == tag:
                    return float(line.split()[0])
            except IndexError:
                print "Cannot find specified input parameter: ", tag
                return 0

    #---------------------------------------------------------------------------#

    def read_control(self):
        """
        Read control.out to set resolutions.
        """

        import numpy as np
    
        data = np.loadtxt(self.dirname+'/out.tgyro.control')

        self.n_radial     = int(data[0])
        self.n_fields     = int(data[1])
        self.n_iterations = int(data[2])

    #---------------------------------------------------------------------------#

    def fileparser(self,file,spec_num=''):
        """
        Generic parser for standard TGYRO file format.
        NOTE: spec_num refers to ion species index
        """

        import string
        import numpy as np

        data = open(self.dirname+'/'+file,'r').readlines()

        # Data dimensions
        nr = self.n_radial+2
        nc = len(string.split(data[0]))
        nb = self.n_iterations+1

        numdata = np.zeros((nc,nb,nr-2))
        for ib in range(nb):
            try:
                tags=string.split(data[ib*nr])
                null=string.split(data[ib*nr+1])
            except:
                print "WARNING: (data.py) "+file+" shorter than expected."
                return 0

            for ir in range(nr-2):
                row=string.split(data[ib*nr+ir+2])
                for ic in range(nc):
                    numdata[ic,ib,ir] = row[ic]

        # Append ion tags with ion species index
        for ic in range(nc):
            if tags[ic]=='r/a' and spec_num!='':
                continue
            self.data[tags[ic]+str(spec_num)] = numdata[ic,:,:]

    #---------------------------------------------------------------------------#

    def get_residual(self):
        """Read out.tgyro.residual
        """
        import string
        import numpy as np

        fn = 'out.tgyro.residual'
        data = open(self.dirname+'/'+fn,'r').readlines()
        
        # Data dimensions 
        nr = self.n_radial
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
        
