class GYROData:
    """A class for management of GYRO output data.

    Data:

    directory_name = ""
    profile  = {}
    geometry = {}
    t        = {}
    freq     = {}
    diff     = []
    diff_i   = []
    diff_n   = []
    gbflux   = []
    gbflux_i = []
    gbflux_n = []
    moment_u = []
    moment_n = []
    moment_e = []
    moment_v = []
    moment_zero    = []
    flux_velocity  = []
    k_perp_squared = []

    Example Usage:
        >>> from pyrats.gyro.data import GYROData
        >>> import matplotlib.pyplot as plt
        >>> sim = GYROData('example_directory')
        >>> sim.make_gbflux()
        >>> plt.plot(sim.gbflux()[0][0][1])
        >>> plt.show()
    """

    # Methods
    def __init__(self, sim_directory):
        """Constructor reads in data from sim_directory and creates new object.
        """

        self.init_data()
        self.set_directory(sim_directory)
        self.read_profile()
        self.read_geometry()
        self.read_t()

        # The rest of the possible read routines should NOT be done automatically.

        #self.read_freq()
        #self.read_balloon()
        #self.read_gbflux_i()
        #self.read_gbflux_n()
        #self.read_moment_u()
        #self.read_moment_n()
        #self.read_moment_e()
        #self.read_moment_v()
        #self.read_moment_zero()
        #self.read_flux_velocity()
        #self.read_k_perp_squared()

        # There are also routines to make "optional fluxes"

        #self.make_gbflux()
        #self.make_diff()
        #self.make_diff_i()

        # And a routine to equalize the time-length.

        #self.equil_time()

        import sys
        from os.path import expanduser, expandvars
        path = '$GACODE_ROOT/shared/python/pyrats'
        sys.path.append(expanduser(expandvars(path)))

    def __dir__(self):
        return self.loaded


    def init_data(self):
        """Initialize object data."""

        self.directory_name = ""

        self.profile  = {}
        self.geometry = {}
        self.t        = {}
        self.freq     = {}
        self.balloon  = {}

        self.diff           = []
        self.diff_i         = []
        self.diff_n         = []
        self.gbflux         = []
        self.gbflux_i       = []
        self.gbflux_n       = []
        self.moment_u       = []
        self.moment_n       = []
        self.moment_e       = []
        self.moment_v       = []
        self.moment_zero    = []
        self.flux_velocity  = []
        self.k_perp_squared = []
        self.loaded         = []


    def set_directory(self, path):
        """Set the simulation directory."""

        from os.path import expanduser, expandvars
        self.directory_name = expanduser(expandvars(path))


    def get_input(self, input_name):
        """Return the specified variable from input.gyro.gen.

        input_name  -  requested input

        Ex: get_input("TIME_STEP")
        """

        input_file = file(self.directory_name + 'input.gyro.gen', 'r')
        for line in input_file:
            try:
                if line.split()[1] == input_name:
                    return float(line.split()[0])
            except IndexError:
                print "Cannot find specified input parameter: ", input_name
                return 0


    def read_file(self, fname, dSize):
        """Reads GYRO data file.  Returns a numpy array with dimensions
        matching input data."""

        import numpy as np
        import pyrats.fileupload as fileupload 
        import sys
        import os

        fname = 'out.gyro.' + fname
        filename = self.directory_name + '/' + fname
        f = np.array([])
        if fname in os.listdir(self.directory_name):
            f = fileupload.loadtxt(filename, dSize)
            return f
        else:
            return []


    def read_t(self):
        """Read out.gyro.t to get time data."""

        import numpy as np

        try:
            t = np.loadtxt(self.directory_name + '/out.gyro.t')
        except:
            print "ERROR (GYROData): Fatal error!  Missing out.gyro.t."
            sys.exit()

        self.t['n_time']   = len(t[:,0])
        self.t['i_time']   = t[:,0]
        self.t['(c_s/a)t'] = t[:,1]
        self.loaded.append('t')


    def read_freq(self):
        """Reads in frequency data.  Output is dictionary of numpy arrays with
        dimensions: n_n x n_time"""

        import numpy as np
        import sys

        try:
            freq = np.loadtxt(self.directory_name+'/out.gyro.freq').transpose()
        except:
            print "ERROR (GYROData): Missing out.gyro.freq."
            sys.exit()

        temp = freq.reshape((4,self.profile['n_n'],self.t['n_time']), order='F')
     
        self.freq['(a/c_s)w']        = temp[0,:,:]
        self.freq['(a/c_s)gamma']    = temp[1,:,:]
        self.freq['err(a/c_s)w']     = temp[2,:,:]
        self.freq['err(a/c_s)gamma'] = temp[3,:,:]

        self.loaded.append('freq')


    def read_profile(self):
        """Read out.gyro.profile to get control data.  Output is dictionary
        containing necessary information."""

        import sys
        import numpy as np

        try:
            profile = np.loadtxt(self.directory_name + '/out.gyro.profile')
        except:
            print "ERROR (GYROData): Fatal error!  Missing out.gyro.profile."
            sys.exit()

        n_x    = int(profile[0]) 
        n_spec = int(profile[15])

        self.profile['n_x']             = int(profile[0])
        self.profile['n_theta_section'] = int(profile[1])
        self.profile['n_pass']          = int(profile[2])
        self.profile['n_trap']          = int(profile[3])
        self.profile['n_lambda']        = int(profile[2]+profile[3])
        self.profile['n_energy']        = int(profile[4])
        self.profile['n_theta_plot']    = int(profile[5])
        self.profile['n0']              = int(profile[6])
        self.profile['n_n']             = int(profile[7])
        self.profile['d_n']             = int(profile[8])
        self.profile['n_explicit_damp'] = int(profile[9])
        self.profile['nonlinear_flag']  = int(profile[10])
        self.profile['electron_method'] = int(profile[11])
        self.profile['n_field']         = int(profile[12])
        self.profile['n_ion']           = int(profile[13])
        self.profile['n_kinetic']       = int(profile[14])
        self.profile['n_spec']          = int(profile[15])
        self.profile['field_r0_flag']   = int(profile[16])
        self.profile['field_r0_grid']   = int(profile[17])
        self.profile['n_grid_exp']      = int(profile[18])
        self.profile['boundary_method'] = int(profile[19])
        self.profile['r']               = profile[20:(20+n_x)]
        self.profile['q']               = profile[(20+n_x):(20+2*n_x)]
        self.profile['r_s']             = profile[(20+2*n_x):(20+3*n_x)]
        self.profile['q_s']             = profile[(20+3*n_x):(20+4*n_x)]
        # The parameter "mark" is used to keep track of where in the file the
        # program is so that the indicies don't get to convoluted.
        mark = 20 + 4*n_x
        temp = profile[mark:(mark+n_spec*n_x)]
        self.profile['dlntdr_s'] = temp.reshape((n_spec,n_x),order='F')
        mark = mark + n_spec*n_x
        temp = profile[mark:(mark+n_spec*n_x)]
        self.profile['dlnndr_s'] = temp.reshape((n_spec,n_x),order='F')
        mark = mark + n_spec*n_x
        temp = profile[mark:(mark+n_spec*n_x)]
        self.profile['tem_s'] = temp.reshape((n_spec,n_x),order='F')
        mark = mark + n_spec*n_x
        temp = profile[mark:(mark+n_spec*n_x)]
        self.profile['den_s'] = temp.reshape((n_spec,n_x),order='F')
        mark = mark + n_spec*n_x
        self.profile['rmaj_s/r_s'] = profile[mark:(mark+n_x)]
        self.profile['delta_s']    = profile[(mark+n_x):(mark+2*n_x)]
        self.profile['zeta_s']     = profile[(mark+2*n_x):(mark+3*n_x)]
        self.profile['kappa_s']    = profile[(mark+3*n_x):(mark+4*n_x)]
        self.profile['drmaj_s']    = profile[(mark+4*n_x):(mark+5*n_x)]
        self.profile['shat_s']     = profile[(mark+5*n_x):(mark+6*n_x)]
        self.profile['s_delta_s']  = profile[(mark+6*n_x):(mark+7*n_x)]
        self.profile['s_zeta_z']   = profile[(mark+7*n_x):(mark+8*n_x)]
        self.profile['s_kappa_s']  = profile[(mark+8*n_x):(mark+9*n_x)]
        self.profile['zmag_s']     = profile[(mark+9*n_x):(mark+10*n_x)]
        self.profile['dzmag_s']    = profile[(mark+10*n_x):(mark+11*n_x)]
        self.profile['beta_unit_s'] = profile[(mark+11*n_x):(mark+12*n_x)]
        self.profile['gamma_e_s']  = profile[(mark+12*n_x):(mark+13*n_x)]
        self.profile['gamma_p_s']  = profile[(mark+13*n_x):(mark+14*n_x)]
        self.profile['mach_s']     = profile[(mark+14*n_x):(mark+15*n_x)]
        self.profile['b_unit_s']   = profile[(mark+15*n_x):(mark+16*n_x)]
        self.profile['dr_eodr']    = profile[(mark+16*n_x):(mark+17*n_x)]
        self.profile['z_eff_s']    = profile[(mark+17*n_x):(mark+18*n_x)]
        self.profile['nu_s']       = profile[(mark+18*n_x):(mark+19*n_x)]
        self.profile['w0_s']       = profile[(mark+19*n_x):(mark+20*n_x)]
        mark = mark + 20*n_x
        self.profile['box_multiplier'] = profile[mark]
        self.profile['lambda']     = profile[(mark+1):(mark+1+self.profile['n_lambda'])]
        mark = mark + 1 + self.profile['n_lambda']
        self.profile['energy']     = profile[mark:(mark+self.profile['n_energy'])]
        self.profile['lambda_tp']  = profile[mark + self.profile['n_energy']]
        mark = mark + self.profile['n_energy'] + 1
        self.profile['kt_rho']     = profile[mark:(mark+self.profile['n_n'])]
        self.profile['rho_s']      = profile[mark+self.profile['n_n']]
        mark = mark + self.profile['n_n'] + 1
        self.profile['z']          = profile[mark:(mark+n_spec)]
        self.profile['n_fine']     = int(profile[mark+n_spec])
        self.profile['n_moment']   = int(profile[mark+n_spec+1])

        if self.profile['n_theta_plot'] == 1:
            self.profile['theta_plot'] = 0.0
        else:
            n_theta = self.profile['n_theta_plot']
            self.profile['theta_plot'] = -np.pi+2*np.pi*np.arange(n_theta)/float(n_theta)


    def read_geometry(self):
        """Reads in geometry_array data.  Output is dictionary of numpy arrays
        with dimensions: n_fine x n_x."""
        
        geometry = self.read_file('geometry_arrays', 16)
        if len(geometry) > 0:
            temp = geometry.reshape((11, self.profile['n_fine'], self.profile['n_x']), order='F')
            self.geometry['v']       = temp[0,:,:]
            self.geometry['gsin']    = temp[1,:,:]
            self.geometry['gcos1']   = temp[2,:,:]
            self.geometry['gcos2']   = temp[3,:,:]
            self.geometry['usin']    = temp[4,:,:]
            self.geometry['ucos']    = temp[5,:,:]
            self.geometry['B']       = temp[6,:,:]
            self.geometry['G_theta'] = temp[7,:,:]
            self.geometry['grad_r']  = temp[8,:,:]
            self.geometry['G_q']     = temp[9,:,:]
            self.geometry['THETA']   = temp[10,:,:]


    def read_gbflux_i(self):
        """Reads in gbflux_i data.  Output is numpy array with dimensions:
        n_kinetic x n_field x 4 x n_x x n_time"""

        import sys

        n_kinetic = self.profile['n_kinetic']
        n_field   = self.profile['n_field']
        n_x       = self.profile['n_x']
      
        gbflux_i = self.read_file('gbflux_i',12)
        if len(gbflux_i) > 0:
            t = len(gbflux_i)/(n_kinetic*n_field*4*n_x)
            self.gbflux_i = gbflux_i.reshape((n_kinetic,n_field,4,n_x,t),order='F')
            self.loaded.append('gbflux_i')
        else:
            print "ERROR (GYROData): out.gyro.gbflux_i not found."
            sys.exit()


    def read_gbflux_n(self):
        """Reads gbflux_n data.  Output is numpy array with dimensions:
        n_kinetic x n_field x 4 x n_n x n_time"""

        import sys

        n_kinetic = self.profile['n_kinetic']
        n_field   = self.profile['n_field']
        n_n       = self.profile['n_n']

        gbflux_n = self.read_file('gbflux_n',12)
        if len(gbflux_n) > 0:
            t = len(gbflux_n)/(n_kinetic*n_field*4*n_n)
            self.gbflux_n = gbflux_n.reshape((n_kinetic,n_field,4,n_n,t),order='F')
            self.loaded.append('gbflux_n')
        else:
            print "ERROR (GYROData): out.gyro.gbflux_n not found."
            sys.exit()


    def read_moment_u(self):
        """Reads in moment_u data.  Output is numpy array with dimensions:
        n_theta_plot x n_x x n_field x n_n x n_time"""

        n_theta_plot = self.profile['n_theta_plot']
        n_x          = self.profile['n_x']
        n_field      = self.profile['n_field']
        n_n          = self.profile['n_n']

        moment_u = self.read_file('moment_u',12)
        if len(moment_u) > 0:
            t = len(moment_u)/(2*n_theta_plot*n_x*n_field*n_n)
            self.moment_u = moment_u.reshape((2,n_theta_plot,n_x,n_field,n_n,t),order='F')
            self.moment_u = self.moment_u[0] + 1j*self.moment_u[1]
            self.loaded.append('moment_u')


    def read_moment_n(self):
        """Reads in moment_n data.  Output is numpy array with dimensions:
        n_theta_plot x n_x x n_kinetic x n_n x n_time"""

        n_theta_plot = self.profile['n_theta_plot']
        n_x          = self.profile['n_x']
        n_kinetic    = self.profile['n_kinetic']
        n_n          = self.profile['n_n']

        moment_n = self.read_file('moment_n',12)
        if len(moment_n) > 0:
            t = len(moment_n)/(2*n_theta_plot*n_x*n_kinetic*n_n)
            self.moment_n = moment_n.reshape((2,n_theta_plot,n_x,n_kinetic,n_n,t),order='F')
            self.moment_n = self.moment_n[0] + 1j*self.moment_n[1]
            self.loaded.append('moment_n')


    def read_moment_e(self):
        """Reads in moment_e data.  Output is numpy array with dimensions:
        n_theta_plot x n_x x n_kinetic x n_n x n_time"""

        n_theta_plot = self.profile['n_theta_plot']
        n_x          = self.profile['n_x']
        n_kinetic    = self.profile['n_kinetic']
        n_n          = self.profile['n_n']

        moment_e = self.read_file('moment_e',12)
        if len(moment_e) > 0:
            t = len(moment_e)/(2*n_theta_plot*n_x*n_kinetic*n_n)
            self.moment_e = moment_e.reshape((2,n_theta_plot,n_x,n_kinetic,n_n,t),order='F')
            self.moment_e = self.moment_e[0] + 1j*self.moment_e[1]
            self.loaded.append('moment_e')


    def read_moment_v(self):
        """Reads in moment_v data.  Output is numpy array with dimensions:
        n_theta_plot x n_x x n_kinetic x n_n x n_time"""

        n_theta_plot = self.profile['n_theta_plot']
        n_x          = self.profile['n_x']
        n_kinetic    = self.profile['n_kinetic']
        n_n          = self.profile['n_n']

        moment_v = self.read_file('moment_v',12)
        if len(moment_v) > 0:
            t = len(moment_v)/(2*n_theta_plot*n_x*n_kinetic*n_n)
            self.moment_v = moment_v.reshape((2,n_theta_plot,n_x,n_kinetic,n_n,t),order='F')
            self.moment_v = self.moment_v[0] + 1j*self.moment_v[1]
            self.loaded.append('moment_v')


    def read_moment_zero(self):
        """Reads in moment_zero data.  Output is numpy array with dimensions:
        n_x x n_kinetic x n_moment x n_time"""

        n_x          = self.profile['n_x']
        n_kinetic    = self.profile['n_kinetic']
        n_moment     = self.profile['n_moment']

        moment_zero = self.read_file('moment_zero',12)
        if len(moment_zero) > 0:
            t = len(moment_zero)/(n_x*n_kinetic*n_moment)
            self.moment_zero = moment_zero.reshape((n_x,n_kinetic,n_moment,t),order='F')
            self.loaded.append('moment_zero')


    def read_flux_velocity(self):
        """Reads out.gyro.flux_velocity.  
        Output is numpy array with dimensions: (n_energy,n_lambda,n_kinetic,n_field,2,n_n,n_time)"""
        n_energy  = self.profile['n_energy']
        n_lambda  = self.profile['n_lambda']
        n_kinetic = self.profile['n_kinetic']
        n_field   = self.profile['n_field']
        n_n       = self.profile['n_n']

        flux_velocity = self.read_file('flux_velocity', 12)
        if len(flux_velocity) > 0:
            t = len(flux_velocity)/(n_energy*n_lambda*n_kinetic*n_field*2*n_n)
            self.flux_velocity = flux_velocity.reshape(
                (n_energy,n_lambda,n_kinetic,n_field,2,n_n,t),order='F')
            self.loaded.append('flux_velocity')


    def read_k_perp_squared(self):
        """Reads out.gyro.k_perp_squared.  
           Output is numpy array with dimensions: (n_n,n_time)"""

        k_perp_squared = self.read_file('k_perp_squared', 12)
        if len(k_perp_squared) > 0:
            t = len(k_perp_squared)/self.profile['n_n']
            self.k_perp_squared = k_perp_squared.reshape((self.profile['n_n'],t),order='F')
            self.loaded.append('k_perp_squared')

    def read_balloon(self):
        """Reads out.gyro.balloon*.  Data is stored in self.balloon"""

        import glob
        import string
        import sys

        m     = self.profile['box_multiplier']
        n_x   = self.profile['n_x']
        n_ang = self.profile['n_theta_plot']*n_x/m

        list = glob.glob(self.directory_name+'/out.gyro.balloon*')

        # If list is empty, then exit with error message.
        if len(list) == 0:
            print "ERROR (GYROData): out.gyro.balloon* not found."
            sys.exit()

        for filename in list:
            ext = string.splitfields(filename,'.')[-1]
            data = self.read_file(ext,12)
            u = data.reshape((2,n_ang,m,self.t['n_time']),order='F') 
            self.balloon[ext] = u[0,...]+1j*u[1,...]
            
            

    #------------------------------------------------------------
    # Create data from other previously imported data

    def make_gbflux(self):
        """Makes gbflux.  Output is numpy array with dimensions:
        n_kinetic x n_field x 4 x n_time"""

        import numpy as np

        if self.profile['boundary_method'] == 1:
            # Periodic simulation
            self.gbflux = np.mean(self.gbflux_i, axis=3)
        else:
            # Nonperiodic simulation: don't include buffers in average
            n = self.profile['n_explicit_damp']
            self.gbflux = np.mean(self.gbflux_i[:,:,:,n:-n,:],axis=3)


    def make_diff(self):
        """Makes diff.  Output is dictionary of numpy arrays with
        dimensions: n_kinetic x n_field x n_time"""

        # NOTE: deprecate n_x_offset in GYRO.

        import numpy as np

        if self.gbflux == []:
            self.make_gbflux()

        ir_norm = int(self.profile['n_x']/2+1)

        n_kinetic = self.profile['n_kinetic']
        n_field   = self.profile['n_field']
        n_time    = self.t['n_time']

        self.diff = np.zeros((n_kinetic,n_field,2,n_time))
        for i in range(n_kinetic):
            # Density diffusivity
            self.diff[i,:,0,:] = self.gbflux[i,:,0,:]/(
                self.profile['dlnndr_s'][i,ir_norm]*
                self.profile['den_s'][i,ir_norm])
            # Energy diffusivity
            self.diff[i,:,1,:] = self.gbflux[i,:,1,:]/(
                self.profile['dlntdr_s'][i,ir_norm]*
                self.profile['tem_s'][i,ir_norm]*
                self.profile['den_s'][i,ir_norm])


    def make_diff_i(self):
        """Makes diff_i.  Output is dictionary of numpy arrays with
        dimensions: n_kinetic x n_field x n_x x n_time"""

        import numpy as np

        if self.gbflux_i == []:
            self.read_gbflux_i()
 
        ir_norm = int(self.profile['n_x']/2+1)

        n_kinetic = self.profile['n_kinetic']
        n_field   = self.profile['n_field']
        n_x       = self.profile['n_x']
        n_time    = self.t['n_time']

        self.diff_i = np.zeros((n_kinetic,n_field,2,n_x,n_time))
        for i in range(n_kinetic):
            # Density diffusivity
            self.diff_i[i,:,0,:,:] = self.gbflux_i[i,:,0,:,:]/(
                self.profile['dlnndr_s'][i,ir_norm]*
                self.profile['den_s'][i,ir_norm])
            # Energy diffusivity
            self.diff_i[i,:,1,:,:] = self.gbflux_i[i,:,1,:,:]/(
                self.profile['dlntdr_s'][i,ir_norm]*
                self.profile['tem_s'][i,ir_norm]*
                self.profile['den_s'][i,ir_norm])


    def plot(self, x, y, dim=(1,1)):
        """Creates matplotlib plot of requested data."""

        import matplotlib.pyplot as plt
        import numpy as np

        fig = plt.figure(self.fignum)
        self.fignum = self.fignum + 1
        ax = fig.add_subplot(dim[0], dim[1], self.plotcounter)
        ax.plot(x, y)

    #----------------------------------------------------------------------#

    def equil_time(self):
        """Equalizes the lengths of the different data arrays, in the case that
        the time axis is longer in some than in others."""

        import numpy as np

        temp = []

        for item in self.loaded:
            if item == 't':
                temp.append(self.t['n_time'])
            elif item == 'freq':
                temp.append(len(self.freq['(a/c_s)w'][-1]))
            else:
                temp.append(len(eval('self.' + item).T))

        cutoff = min(temp)

        for item in self.loaded:

            if item == 't':
                n = self.t['n_time']
                self.t['i_time']   = self.t['i_time'][...,0:cutoff]
                self.t['(c_s/a)t'] = self.t['(c_s/a)t'][...,0:cutoff]
                self.t['n_time']   = cutoff
            elif item == 'freq':
                self.freq['(a/c_s)w']        = self.freq['(a/c_s)w'][...,0:cutoff]
                self.freq['(a/c_s)gamma']    = self.freq['(a/c_s)gamma'][...,0:cutoff]
                self.freq['err(a/c_s)w']     = self.freq['err(a/c_s)w'][...,0:cutoff]
                self.freq['err(a/c_s)gamma'] = self.freq['err(a/c_s)gamma'][...,0:cutoff]
            elif item == 'gbflux_i':
                self.gbflux_i = self.gbflux_i[...,0:cutoff] 
            elif item == 'gbflux_n':
                self.gbflux_n = self.gbflux_n[...,0:cutoff]
            elif item == 'moment_u':
                self.moment_u = self.moment_u[...,0:cutoff]
            elif item == 'moment_n':
                self.moment_n = self.moment_n[...,0:cutoff]
            elif item == 'moment_e':
                self.moment_e = self.moment_e[...,0:cutoff]
            elif item == 'moment_v':
                self.moment_v = self.moment_v[...,0:cutoff]
            elif item == 'moment_zero':
                self.moment_zero = self.moment_zero[...,0:cutoff]
            elif item == 'flux_velocity':
                self.flux_velocity = self.flux_velocity[...,0:cutoff]
            elif item == 'k_perp_squared':
                self.k_perp_squared = self.k_perp_squared[...,0:cutoff]
