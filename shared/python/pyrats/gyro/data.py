class GYROData:
    """GYRO output data class.

    Data:
    
    dirname = ""
    profile  = {}
    geometry = {}
    t        = {}
    freq     = {}
    balloon  = []
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
        >>> sim = GYROData('example_directory')
        >>> sim.make_gbflux()
    """

    #---------------------------------------------------------------------------#
    # Methods

    def __init__(self, sim_directory):
        """Constructor reads in data from sim_directory and creates new object.
        """

        self.init_data()
        self.set_directory(sim_directory)
        self.read_profile()
        self.read_geometry()
        self.read_t()
        self.make_tags()
        self.debug=False

        # The rest of the possible read routines should NOT be done automatically.

        #self.read_freq()
        #self.read_balloon()
        #self.read_gbflux_i()
        #self.read_gbflux_n()
        #self.read_gbflux_exc()
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

        import sys, os
        from os.path import expanduser, expandvars
        path = '$GACODE_ROOT/shared/python/pyrats'
        sys.path.append(expanduser(expandvars(path)))

    #---------------------------------------------------------------------------#

    def init_data(self):
        """Initialize object data."""

        self.dirname = ""
        self.debug=False

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
        self.gbflux_exc     = []
        self.moment_u       = []
        self.moment_n       = []
        self.moment_e       = []
        self.moment_v       = []
        self.moment_zero    = []
        self.flux_velocity  = []
        self.k_perp_squared = []
        self.loaded         = []

        self.tagspec    = []
        self.tagmom     = []
        self.tagmomtext = []
 
    #---------------------------------------------------------------------------#
        
    def set_directory(self, path):
        """Set the simulation directory."""

        from os.path import expanduser, expandvars
        self.dirname = expanduser(expandvars(path))

    #---------------------------------------------------------------------------#

    def read_t(self):
        """Read out.gyro.t to get time data."""

        import sys, os
        import numpy as np

        try:
            import tables

            h5file=self.dirname+'/out.gyro.timedata.h5'
            ascfile=self.dirname+'/out.gyro.t'
            if os.path.exists(h5file):
              f=tables.openFile(h5file,'r')
              self.t['data_step']   = f.root.data_step.read()
              self.t['(c_s/a)t'] = f.root.t_current.read()
              self.t['n_time']   = self.t['i_time'].shape[0]
              f.close()
            else:
              print "ERROR (GYROData): File does not exist "+h5file
              sys.exit()
        except:
            if self.debug:
                print "Unable to use HDF5 data"
         
        print os.path.exists(ascfile)
        if os.path.exists(ascfile):
            try:
                t = np.loadtxt(ascfile)
                self.t['n_time']   = len(t[:,0])
                self.t['data_step']   = t[:,0]
                self.t['(c_s/a)t'] = t[:,1]
            except:
                print "ERROR (GYROData): Problem loading data from "+ascfile
                sys.exit()
        else:
          print "ERROR (GYROData): Fatal error!  Missing "+ascfile
          sys.exit()

        self.loaded.append('t')

        self.n = len(t[:,0])

    #---------------------------------------------------------------------------#

    def read_freq(self):
        """Reads in frequency data.  Output is dictionary of numpy arrays with
        dimensions: n_n x n_time"""

        import numpy as np
        import sys

        h5file=self.dirname+'/out.gyro.timedata.h5'
        ascfile=self.dirname+'/out.gyro.freq'
        if os.path.exists(h5file):
            try:
                print "Frequency data not implemented in HDF5 yet"
                f=tables.openFile(h5file,'r')
                self.freq['(a/c_s)w']        = f.root.omega.read()
                self.freq['(a/c_s)gamma']    = f.root.gamma.read()
                self.freq['err(a/c_s)w']     = f.root.error_omega.read()
                self.freq['err(a/c_s)gamma'] = f.root.error_gamma.read()
                f.close()
            except:
                print "ERROR (GYROData): Problem reading "+h5file
                sys.exit()
        if os.path.exists(ascfile):
            try:
                freq = np.loadtxt(ascfile).transpose()
                self.freq['(a/c_s)w']        = temp[0,:,:]
                self.freq['(a/c_s)gamma']    = temp[1,:,:]
                self.freq['err(a/c_s)w']     = temp[2,:,:]
                self.freq['err(a/c_s)gamma'] = temp[3,:,:]
                temp = freq.reshape((4,self.profile['n_n'],self.t['n_time']), order='F')
            except:
              if self.debug:
                  print "Unable to use HDF5 data"
        else:
            print "ERROR (GYROData): Missing out.gyro.freq."
            sys.exit()


        self.loaded.append('freq')

    #---------------------------------------------------------------------------#

    def read_profile(self):
        """Read out.gyro.profile to get control data.  Output is dictionary
        containing necessary information."""

        import sys, os
        import numpy as np

        h5file=self.dirname+'/out.gyro.initdata.h5'
        ascfile=self.dirname+'/out.gyro.profile'
        try:
          import tables

          if os.path.exists(h5file):
              try:
                  f=tables.openFile(h5file,'r')
                  self.profile['n_x']             = f.root.n_x.read()
                  self.profile['n_theta_section'] = f.root.n_theta_section.read()
                  self.profile['n_pass']          = f.root.n_pass.read()
                  self.profile['n_trap']          = f.root.n_trap.read() 
                  self.profile['n_lambda']        = self.profile['n_pass']+self.profile['n_trap']
                  self.profile['n_energy']        = f.root.n_energy.read() 
                  self.profile['n_theta_plot']    = f.root.n_theta_plot.read() 
                  self.profile['n0']              = f.root.n0.read() 
                  self.profile['n_n']             = f.root.n_n.read() 
                  self.profile['d_n']             = f.root.d_n.read() 
                  self.profile['n_explicit_damp'] = f.root.n_explicit_damp.read() 
                  self.profile['nonlinear_flag']  = f.root.nonlinear_flag.read() 
                  self.profile['electron_method'] = f.root.electron_method.read() 
                  self.profile['n_field']         = f.root.n_field.read() 
                  self.profile['n_ion']           = f.root.n_ion.read() 
                  self.profile['n_kinetic']       = f.root.n_kinetic.read() 
                  self.profile['n_spec']          = f.root.n_spec.read() 
                  self.profile['field_r0_flag']   = f.root.field_r0_flag.read() 
                  self.profile['field_r0_grid']   = f.root.field_r0_grid.read() 
                  self.profile['n_grid_exp']      = f.root.n_grid_exp.read() 
                  self.profile['boundary_method'] = f.root.boundary_method.read() 
                  self.profile['r']               = f.root.r.read() 
                  self.profile['q']               = f.root.q.read() 
                  self.profile['r_s']             = f.root.r_s.read() 
                  self.profile['q_s']             = f.root.q_s.read() 
                  self.profile['dlntdr_s']        = f.root.dlntdr_s.read() 
                  self.profile['dlnndr_s']        = f.root.dlnndr_s.read() 
                  self.profile['tem_s']           = f.root.tem_s.read() 
                  self.profile['den_s']           = f.root.den_s.read() 
                  self.profile['rmaj_s/r_s']      = f.root.rmaj_s.read() 
                  self.profile['delta_s']         = f.root.delta_s.read() 
                  self.profile['zeta_s']          = f.root.zeta_s.read() 
                  self.profile['kappa_s']         = f.root.kappa_s.read() 
                  self.profile['drmaj_s']         = f.root.drmaj_s.read() 
                  self.profile['shat_s']          = f.root.shat_s.read() 
                  self.profile['s_delta_s']       = f.root.s_delta_s.read() 
                  self.profile['s_zeta_s']        = f.root.s_zeta_s.read() 
                  self.profile['s_kappa_s']       = f.root.s_kappa_s.read() 
                  self.profile['zmag_s']          = f.root.zmag_s.read() 
                  self.profile['dzmag_s']         = f.root.dzmag_s.read() 
                  self.profile['beta_unit_s']     = f.root.beta_unit_s.read() 
                  self.profile['gamma_e_s']       = f.root.gamma_e_s.read() 
                  self.profile['gamma_p_s']       = f.root.gamma_p_s.read() 
                  self.profile['mach_s']          = f.root.mach_s.read() 
                  self.profile['b_unit_s']        = f.root.b_unit_s.read() 
                  self.profile['dr_eodr']         = f.root.dr_eodr.read() 
                  self.profile['z_eff_s']         = f.root.z_eff_s.read() 
                  self.profile['nu_s']            = f.root.nu_s.read() 
                  self.profile['w0_s']            = f.root.w0_s.read() 
                  self.profile['box_multiplier']  = f.root.box_multiplier.read() 
                  #self.profile['lambda']          = eval("f.root.lambda.read()") 
                  self.profile['energy']          = f.root.energy.read() 
                  self.profile['lambda_tp']       = f.root.lambda_tp.read() 
                  self.profile['kt_rho']          = f.root.kt_rho.read() 
                  self.profile['rho_s']           = f.root.rho_s.read() 
                  self.profile['z']               = f.root.z.read() 
                  self.profile['n_fine']          = f.root.n_fine.read() 
                  self.profile['n_moment']        = f.root.n_moment.read() 
                  if self.profile['n_theta_plot'] == 1:
                      self.profile['theta_plot'] = 0.0
                  else:
                      n_theta = self.profile['n_theta_plot']
                      self.profile['theta_plot'] = -np.pi+2*np.pi*np.arange(n_theta)/float(n_theta)
                  f.close()

              except:
                  print "ERROR (GYROData): Fatal error! Missing "+h5file
                  raise
                  sys.exit()
        except:
            print "Problem importing h5 file"
            raise

        if os.path.exists(ascfile):
            try:
                profile = np.loadtxt(ascfile)
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
                self.profile['s_zeta_s']   = profile[(mark+7*n_x):(mark+8*n_x)]
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

            except:
                print "ERROR (GYROData): Fatal error! While reading "+ascfile
                sys.exit()
        else:
            print "ERROR (GYROData): Fatal error!  Missing "+ascfile
            sys.exit()


    #---------------------------------------------------------------------------#

    def read_geometry(self):
        """Reads in geometry_array data.  Output is dictionary of numpy arrays
        with dimensions: n_fine x n_x."""
        
        import sys, os
        import numpy as np

        h5file=self.dirname+'/out.gyro.initdata.h5'
        ascfile=self.dirname+'/out.gyro.t'
        if os.path.exists(h5file):
            try:
                f=tables.openFile(h5file,'r')
                self.geometry['v']       = f.root.v.read() 
                self.geometry['gsin']    = f.root.gsin.read() 
                self.geometry['gcos1']   = f.root.gcos1.read() 
                self.geometry['gcos2']   = f.root.gcos2.read() 
                self.geometry['usin']    = f.root.usin.read() 
                self.geometry['ucos']    = f.root.ucos.read() 
                self.geometry['B']       = f.root.B.read() 
                self.geometry['G_theta'] = f.root.G_theta.read() 
                self.geometry['grad_r']  = f.root.grad_r.read() 
                self.geometry['G_q']     = f.root.G_q.read() 
                self.geometry['THETA']   = f.root.THETA.read() 
                f.close()
            except:
                print "ERROR (GYROData): File does not exist "+h5file
                sys.exit()
        if os.path.exists(ascfile):
            try:
                geometry = np.fromfile(self.dirname+'/out.gyro.geometry_arrays',dtype=float,sep=" ")
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
            except:
                print "ERROR (GYROData): File does not exist "+ascfile
                sys.exit()
        else:
            print "ERROR (GYROData): Fatal error!  Missing "+ascfile
            sys.exit()

        self.loaded.append('geometry')

    #---------------------------------------------------------------------------#

    def read_gbflux_i(self):
        """Reads in gbflux_i data.  Output is numpy array with dimensions:
        n_kinetic x n_field x 4 x n_x x n_time"""

        import sys, os
        import numpy as np

        n_kinetic = self.profile['n_kinetic']
        n_field   = self.profile['n_field']
        n_x       = self.profile['n_x']
      
        h5file=self.dirname+'/out.gyro.timedata.h5'
        ascfile=self.dirname+'/out.gyro.gbflux_i'
        if os.path.exists(h5file):
            try:
                f=tables.openFile(h5file,'r')
                self.gbflux_i = f.root.gbflux_i.read()
                f.close()
            except:
                print "ERROR (GYROData): File does not exist "+h5file
                sys.exit()
        if os.path.exists(ascfile):
            try:
                gbflux_i = np.fromfile(ascfile,dtype=float,sep=" ")
                nt = len(gbflux_i)/(n_kinetic*n_field*4*n_x)

                if self.n > nt:
                    print 'ERROR (GYROData): '+self.dirname+'/out.gyro.gbflux_i too small. '
                    sys.exit()
                
                self.gbflux_i = gbflux_i.reshape((n_kinetic,n_field,4,n_x,nt),order='F')
            except:
                print "ERROR (GYROData): File does not exist "+ascfile
                sys.exit()
        else:
            print "ERROR (GYROData): Fatal error!  Missing "+ascfile
            sys.exit()

        self.loaded.append('gbflux_i')

    #---------------------------------------------------------------------------#

    def read_gbflux_n(self):
        """Reads gbflux_n data.  Output is numpy array with dimensions:
        n_kinetic x n_field x 4 x n_n x n_time"""

        import sys, os
        import numpy as np

        n_kinetic = self.profile['n_kinetic']
        n_field   = self.profile['n_field']
        n_n       = self.profile['n_n']

        h5file=self.dirname+'/out.gyro.timedata.h5'
        ascfile=self.dirname+'/out.gyro.gbflux_n'
        if os.path.exists(h5file):
            try:
                f=tables.openFile(h5file,'r')
                self.gbflux_i = f.root.gbflux_n.read()
                f.close()
            except:
                print "ERROR (GYROData): File does not exist "+h5file
                sys.exit()
        if os.path.exists(ascfile):
            try:
                gbflux_n = np.fromfile(self.dirname+'/out.gyro.gbflux_n',dtype=float,sep=" ")
                nt = len(gbflux_n)/(n_kinetic*n_field*4*n_n)

                if self.n > nt:
                    print 'ERROR (GYROData): '+self.dirname+'/out.gyro.gbflux_n too small. '
                    sys.exit()
                
                self.gbflux_n = gbflux_n.reshape((n_kinetic,n_field,4,n_n,nt),order='F')
            except:
                print "ERROR (GYROData): File does not exist "+ascfile
                sys.exit()
        else:
            print "ERROR (GYROData): Fatal error!  Missing "+ascfile
            sys.exit()

        self.loaded.append('gbflux_n')

     #---------------------------------------------------------------------------#

    def read_gbflux_exc(self):
        """Reads gbflux_exc data.  Output is numpy array with dimensions:
        n_kinetic x 4 x n_time"""

        import sys, os
        import numpy as np

        n_kinetic = self.profile['n_kinetic']

        try:
            gbflux_exc = np.fromfile(self.dirname+'/out.gyro.gbflux_exc',dtype=float,sep=" ")
        except:
            print "ERROR (GYROData): out.gyro.gbflux_exc not found."
            sys.exit()

        nt = len(gbflux_exc)/(n_kinetic*4)

        if self.n > nt:
            print 'ERROR (GYROData): '+self.dirname+'/out.gyro.gbflux_exc too small. '
            sys.exit()
        
        self.gbflux_exc = gbflux_exc.reshape((n_kinetic,4,nt),order='F')
        self.loaded.append('gbflux_exc')

    #---------------------------------------------------------------------------#

    def read_moment_u(self):
        """Reads in moment_u data.  Output is numpy array with dimensions:
        n_theta_plot x n_x x n_field x n_n x n_time"""

        import sys, os
        import numpy as np

        n_theta_plot = self.profile['n_theta_plot']
        n_x          = self.profile['n_x']
        n_field      = self.profile['n_field']
        n_n          = self.profile['n_n']

        h5file=self.dirname+'/out.gyro.timedata.h5'
        ascfile=self.dirname+'/out.gyro.moment_u'
        if os.path.exists(h5file):
            try:
                f=tables.openFile(h5file,'r')
                self.moment_u = f.root.moment_uread()
                f.close()
            except:
                print "ERROR (GYROData): File does not exist "+h5file
                sys.exit()
        if os.path.exists(ascfile):
            try:
                data = np.fromfile(ascfile,dtype=float,sep=" ")
                nt = len(data)/(2*n_theta_plot*n_x*n_field*n_n)

                if self.n > nt:
                    print 'ERROR (GYROData): '+self.dirname+'/out.gyro.moment_u too small. '
                    sys.exit()

                self.moment_u = data.reshape((2,n_theta_plot,n_x,n_field,n_n,nt),order='F')
                self.moment_u = self.moment_u[0] + 1j*self.moment_u[1]

            except:
                print "ERROR (GYROData): File does not exist "+ascfile
                sys.exit()
        else:
            print "ERROR (GYROData): Fatal error!  Missing "+ascfile
            sys.exit()

        self.loaded.append('moment_u')


    #---------------------------------------------------------------------------#

    def read_moment_n(self):
        """Reads in moment_n data.  Output is numpy array with dimensions:
        n_theta_plot x n_x x n_kinetic x n_n x n_time"""

        import sys, os
        import numpy as np

        n_theta_plot = self.profile['n_theta_plot']
        n_x          = self.profile['n_x']
        n_kinetic    = self.profile['n_kinetic']
        n_n          = self.profile['n_n']

        h5file=self.dirname+'/out.gyro.timedata.h5'
        ascfile=self.dirname+'/out.gyro.moment_n'
        if os.path.exists(h5file):
            try:
                f=tables.openFile(h5file,'r')
                self.moment_n = f.root.moment_n.read()
                f.close()
            except:
                print "ERROR (GYROData): File does not exist "+h5file
                sys.exit()
        if os.path.exists(ascfile):
            try:
                data = np.fromfile(ascfile,dtype=float,sep=" ")

                nt = len(data)/(2*n_theta_plot*n_x*n_kinetic*n_n)

                if self.n > nt:
                    print 'ERROR (GYROData): '+self.dirname+'/out.gyro.moment_n too small. '
                    sys.exit()

                self.moment_n = data.reshape((2,n_theta_plot,n_x,n_kinetic,n_n,nt),order='F')
                self.moment_n = self.moment_n[0] + 1j*self.moment_n[1]

            except:
                print "ERROR (GYROData): File does not exist "+ascfile
                sys.exit()
        else:
            print "ERROR (GYROData): Fatal error!  Missing "+ascfile
            sys.exit()

        self.loaded.append('moment_n')

    #---------------------------------------------------------------------------#

    def read_moment_e(self):
        """Reads in moment_e data.  Output is numpy array with dimensions:
        n_theta_plot x n_x x n_kinetic x n_n x n_time"""

        import sys, os
        import numpy as np

        n_theta_plot = self.profile['n_theta_plot']
        n_x          = self.profile['n_x']
        n_kinetic    = self.profile['n_kinetic']
        n_n          = self.profile['n_n']

        h5file=self.dirname+'/out.gyro.timedata.h5'
        ascfile=self.dirname+'/out.gyro.moment_e'
        if os.path.exists(h5file):
            try:
                f=tables.openFile(h5file,'r')
                self.moment_e = f.root.moment_e.read()
                f.close()
            except:
                print "ERROR (GYROData): File does not exist "+h5file
                sys.exit()
        if os.path.exists(ascfile):
            try:
                data = np.fromfile(ascfile,dtype=float,sep=" ")

                nt = len(data)/(2*n_theta_plot*n_x*n_kinetic*n_n)

                if self.n > nt:
                    print 'ERROR (GYROData): '+self.dirname+'/out.gyro.moment_e too small. '
                    sys.exit()

                self.moment_e = data.reshape((2,n_theta_plot,n_x,n_kinetic,n_n,nt),order='F')
                self.moment_e = self.moment_e[0] + 1j*self.moment_e[1]

            except:
                print "ERROR (GYROData): File does not exist "+ascfile
                sys.exit()
        else:
            print "ERROR (GYROData): Fatal error!  Missing "+ascfile
            sys.exit()

        self.loaded.append('moment_e')

    #---------------------------------------------------------------------------#

    def read_moment_v(self):
        """Reads in moment_v data.  Output is numpy array with dimensions:
        n_theta_plot x n_x x n_kinetic x n_n x n_time"""

        import sys, os
        import numpy as np

        n_theta_plot = self.profile['n_theta_plot']
        n_x          = self.profile['n_x']
        n_kinetic    = self.profile['n_kinetic']
        n_n          = self.profile['n_n']

        h5file=self.dirname+'/out.gyro.timedata.h5'
        ascfile=self.dirname+'/out.gyro.moment_v'
        if os.path.exists(h5file):
            try:
                f=tables.openFile(h5file,'r')
                self.moment_v = f.root.moment_v.read()
                f.close()
            except:
                print "ERROR (GYROData): File does not exist "+h5file
                sys.exit()
        if os.path.exists(ascfile):
            try:
                data = np.fromfile(ascfile,dtype=float,sep=" ")
  
                nt = len(data)/(2*n_theta_plot*n_x*n_kinetic*n_n)

                if self.n > nt:
                    print 'ERROR (GYROData): '+self.dirname+'/out.gyro.moment_e too small. '
                    sys.exit()

                self.moment_v = data.reshape((2,n_theta_plot,n_x,n_kinetic,n_n,nt),order='F')
                self.moment_v = self.moment_v[0] + 1j*self.moment_v[1]

            except:
                print "ERROR (GYROData): File does not exist "+ascfile
                sys.exit()
        else:
            print "ERROR (GYROData): Fatal error!  Missing "+ascfile
            sys.exit()

        self.loaded.append('moment_v')

    #---------------------------------------------------------------------------#

    def read_moment_zero(self):
        """Read data in out.gyro.moment_zero, store in self.moment_zero. 
        Dimensions: (n_x,n_kinetic,n_moment,n_time)"""

        import sys, os
        import numpy as np

        n_x          = self.profile['n_x']
        n_kinetic    = self.profile['n_kinetic']
        n_moment     = self.profile['n_moment']

        h5file=self.dirname+'/out.gyro.timedata.h5'
        ascfile=self.dirname+'/out.gyro.gbflux_n'
        if os.path.exists(h5file):
            try:
                f=tables.openFile(h5file,'r')
                self.moment_zero = f.root.moment_zero.read()
                f.close()
            except:
                print "ERROR (GYROData): File does not exist "+h5file
                sys.exit()
        if os.path.exists(ascfile):
            try:
                data = np.fromfile(self.dirname+'/out.gyro.moment_zero',dtype=float,sep=" ")

                t = len(data)/(n_x*n_kinetic*n_moment)
                self.moment_zero = data.reshape((n_x,n_kinetic,n_moment,t),order='F')

            except:
                print "ERROR (GYROData): File does not exist "+ascfile
                sys.exit()
        else:
            print "ERROR (GYROData): Fatal error!  Missing "+ascfile
            sys.exit()

        self.loaded.append('moment_zero')

    #---------------------------------------------------------------------------#

    def read_flux_velocity(self):
        """Reads out.gyro.flux_velocity.  
        Output is numpy array with dimensions: (n_energy,n_lambda,n_kinetic,n_field,2,n_n,n_time)"""
        import sys, os
        import numpy as np

        n_energy  = self.profile['n_energy']
        n_lambda  = self.profile['n_lambda']
        n_kinetic = self.profile['n_kinetic']
        n_field   = self.profile['n_field']
        n_n       = self.profile['n_n']


        h5file=self.dirname+'/out.gyro.timedata.h5'
        ascfile=self.dirname+'/out.gyro.flux_velocity'
        if os.path.exists(h5file):
            try:
                f=tables.openFile(h5file,'r')
                self.flux_velocity = f.root.flux_velocity.read()
                f.close()
            except:
                print "ERROR (GYROData): File does not exist "+h5file
                sys.exit()
        if os.path.exists(ascfile):
            try:
                data = np.fromfile(ascfile,dtype=float,sep=" ")
                t = len(data)/(n_energy*n_lambda*n_kinetic*n_field*2*n_n)
                self.flux_velocity = data.reshape(
                    (n_energy,n_lambda,n_kinetic,n_field,2,n_n,t),order='F')

            except:
                print "ERROR (GYROData): File does not exist "+ascfile
                sys.exit()
        else:
            print "ERROR (GYROData): Fatal error!  Missing "+ascfile
            sys.exit()

        self.loaded.append('flux_velocity')


    #---------------------------------------------------------------------------#

    def read_k_perp_squared(self):
        """Reads out.gyro.k_perp_squared.  
           Output is numpy array with dimensions: (n_n,n_time)"""

        import sys, os
        import numpy as np

        h5file=self.dirname+'/out.gyro.timedata.h5'
        ascfile=self.dirname+'/out.gyro.gbflux_n'
        if os.path.exists(h5file):
            try:
                f=tables.openFile(h5file,'r')
                self.k_perp_squared = f.root.k_perp_squared.read()
                f.close()
            except:
                print "ERROR (GYROData): File does not exist "+h5file
                sys.exit()
        if os.path.exists(ascfile):
            try:
                data = np.fromfile(self.dirname+'/out.gyro.k_perp_squared',dtype=float,sep=" ")
                t = len(data)/self.profile['n_n']
                self.k_perp_squared = data.reshape((self.profile['n_n'],t),order='F')

            except:
                print "ERROR (GYROData): File does not exist "+ascfile
                sys.exit()
        else:
            print "ERROR (GYROData): Fatal error!  Missing "+ascfile
            sys.exit()

        self.loaded.append('k_perp_squared')

    #---------------------------------------------------------------------------#

    def read_balloon(self):
        """Reads out.gyro.balloon*.  Data is stored in self.balloon"""

        import glob
        import string
        import sys, os
        import numpy as np

        m     = self.profile['box_multiplier']
        n_x   = self.profile['n_x']
        n_ang = self.profile['n_theta_plot']*n_x/m

        h5file=self.dirname+'/out.gyro.timedata.h5'
        if os.path.exists(h5file):
            try:
                f=tables.openFile(h5file,'r')
                #INSERT
                f.close()
            except:
                print "ERROR (GYROData): Error reading file "+h5file
                sys.exit()
        else:
            list = glob.glob(self.dirname+'/out.gyro.balloon*')

            # If list is empty, then exit with error message.
            if len(list) == 0:
                print "ERROR (GYROData): out.gyro.balloon* not found."
                sys.exit()

            for filename in list:
                data = np.fromfile(filename,dtype=float,sep=" ")
                u = data.reshape((2,n_ang,m,self.t['n_time']),order='F')
                ext = string.splitfields(filename,'.')[-1]
                self.balloon[ext] = u[0,...]+1j*u[1,...]
            
    #------------------------------------------------------------
    # Create data from other previously imported data

    def make_gbflux(self):
        """Makes gbflux (omitting buffers properly). Output is numpy array
           with dimensions: n_kinetic x n_field x 4 x n_time"""

        import numpy as np

        if self.profile['boundary_method'] == 1:
            # Periodic simulation
            self.gbflux = np.mean(self.gbflux_i, axis=3)
        else:
            # Nonperiodic simulation: don't include buffers in average
            n = self.profile['n_explicit_damp']
            self.gbflux = np.mean(self.gbflux_i[:,:,:,n:-n,:],axis=3)

    #---------------------------------------------------------------------------#

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


    #---------------------------------------------------------------------------#

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


    #---------------------------------------------------------------------------#

    def make_tags(self):
        """Generate tags for fields, moments and species"""                 

        import numpy as np

        n_kinetic = self.profile['n_kinetic']        
        melec     = self.profile['electron_method']

        # Species tags        

        for i in range(n_kinetic):

            if melec == 2 or melec == 4:
                if i == n_kinetic-1:
                    self.tagspec.append('elec')
                else:
                    self.tagspec.append('ion-'+str(i+1))
 
            if melec == 1:
                self.tagspec.append('ion-'+str(i+1))
                    
            if melec == 3:
                self.tagspec.append('elec')
        
       # Moment tags

        self.tagmomtext = ['GAMMA [GB]','Q [GB]','PI [GB]','S [GB]']

        self.tagmom = ['\Gamma/\Gamma_\mathrm{GB}',
                       'Q/Q_\mathrm{GB}',
                       '\Pi/\Pi_\mathrm{GB}',
                       'S/S_\mathrm{GB}']

       # Field tags
        self.tagfieldtext = ['Phi',
                             'Apar',
                             'Bpar',
                             'Tot']        

        self.tagfield = ['\mathrm{electrostatic}',
                         '\mathrm{flutter}',
                         '\mathrm{compression}',
                         '\mathrm{total}']
