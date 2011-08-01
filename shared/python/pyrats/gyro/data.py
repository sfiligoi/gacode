class GYROData:
    """A class of GYRO output data.

    Data:

    directory_name = ""
    profile = {}
    geometry = {}
    t = {}
    freq = {}
    diff = []
    diff_i = []
    diff_n = []
    gbflux = []
    gbflux_i = []
    gbflux_n = []
    moment_u = []
    moment_n = []
    moment_e = []
    moment_v = []
    moment_zero = []
    flux_velocity = []
    k_perp_squared = []

    Example Usage:
        >>> from pyrats.data import GYROData
        >>> import matplotlib.pyplot as plt
        >>> sim1 = GYROData('example_directory')
        >>> sim1.make_gbflux()
        >>> plt.plot(sim1.gbflux()[0][0][1])
        >>> plt.show()
    """

    # Methods
    def __init__(self, sim_directory):
        """Constructor reads in data from sim_directory and creates new object.
        """

        self.init_data()
        self.directory_name = sim_directory
        self.read_data()
        import sys
        from os.path import expanduser, expandvars
        path = '$GACODE_ROOT/shared/python/pyrats'
        sys.path.append(expanduser(expandvars(path)))

    def init_data(self):
        """Initialize object data."""

        self.directory_name = ""
        self.profile = {}
        self.geometry = {}
        self.t = {}
        self.freq = {}
        self.diff = []
        self.diff_i = []
        self.diff_n = []
        self.gbflux = []
        self.gbflux_i = []
        self.gbflux_n = []
        self.moment_u = []
        self.moment_n = []
        self.moment_e = []
        self.moment_v = []
        self.moment_zero = []
        self.flux_velocity = []
        self.k_perp_squared = []
        self.loaded = []

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

    def read_data(self):
        """Read in object data.  Executes read_profile, read_geometry, read_t,
        and read_gbflux_i by default."""
        
        self.read_profile()
        self.read_geometry()
        self.read_t()
        #self.read_freq()
        self.read_gbflux_i()
        self.read_gbflux_n()
        #self.read_moment_u()
        #self.read_moment_n()
        #self.read_moment_e()
        #self.read_moment_v()
        #self.read_moment_zero()
        #self.read_flux_velocity()
        #self.read_k_perp_squared()

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
            print "ERROR: File " + fname + " does not exist."
            return f


    def read_t(self):
        """Read t.out to get time data."""

        import numpy as np

        t = np.loadtxt(file(self.directory_name + '/t.out'))
        self.t['t/deltat'] = t[:, 0]
        self.t['(cbar_s/a)t'] = t[:, 1]
        self.t['n_time'] = len(t[:, 0])

    def read_profile(self):
        """Read out.gyro.profile to get control data.  Output is dictionary
        containing necessary information."""

        import numpy as np
        import math

        profile = np.loadtxt(self.directory_name + '/out.gyro.profile')
        self.profile['n_x'] = profile[0]
        self.profile['n_theta_section'] = profile[1]
        self.profile['n_pass'] = profile[2]
        self.profile['n_trap'] = profile[3]
        self.profile['n_lambda'] = self.profile['n_pass'] + self.profile['n_trap']
        self.profile['n_energy'] = profile[4]
        self.profile['n_theta_plot'] = profile[5]
        self.profile['n0'] = profile[6]
        self.profile['n_n'] = profile[7]
        self.profile['d_n'] = profile[8]
        self.profile['n_explicit_damp'] = profile[9]
        self.profile['nonlinear_flag'] = profile[10]
        self.profile['electron_method'] = profile[11]
        self.profile['n_field'] = profile[12]
        self.profile['n_ion'] = profile[13]
        self.profile['n_kinetic'] = profile[14]
        self.profile['n_spec'] = profile[15]
        self.profile['field_r0_flag'] = profile[16]
        self.profile['field_r0_grid'] = profile[17]
        self.profile['n_grid_exp'] = profile[18]
        self.profile['boundary_method'] = profile[19]
        self.profile['r'] = profile[20:(20 + self.profile['n_x'])]
        self.profile['q'] = profile[(20 + self.profile['n_x']):(20 + 2*self.profile['n_x'])]
        self.profile['r_s'] = profile[(20 + 2*self.profile['n_x']):(20 + 3*self.profile['n_x'])]
        self.profile['q_s'] = profile[(20 + 3*self.profile['n_x']):(20 + 4*self.profile['n_x'])]
        # The parameter "mark" is used to keep track of where in the file the
        # program is so that the indicies don't get to convoluted.
        mark = 20 + 4*self.profile['n_x']
        temp = profile[mark:(mark + self.profile['n_spec']*self.profile['n_x'])]
        self.profile['dlntdr_s'] = temp.reshape( (self.profile['n_spec'], self.profile['n_x']), order='F')
        mark = mark + self.profile['n_spec']*self.profile['n_x']
        temp = profile[mark:(mark + self.profile['n_spec']*self.profile['n_x'])]
        self.profile['dlnndr_s'] = temp.reshape( (self.profile['n_spec'], self.profile['n_x']), order='F')
        mark = mark + self.profile['n_spec']*self.profile['n_x']
        temp = profile[mark:(mark + self.profile['n_spec']*self.profile['n_x'])]
        self.profile['tem_s'] = temp.reshape( (self.profile['n_spec'], self.profile['n_x']), order='F')
        mark = mark + self.profile['n_spec']*self.profile['n_x']
        temp = profile[mark:(mark + self.profile['n_spec']*self.profile['n_x'])]
        self.profile['den_s'] = temp.reshape( (self.profile['n_spec'], self.profile['n_x']), order='F')
        mark = mark + self.profile['n_spec']*self.profile['n_x']
        self.profile['rmaj_s/r_s'] = profile[mark:(mark + self.profile['n_x'])]
        self.profile['delta_s'] = profile[(mark + self.profile['n_x']):(mark + 2*self.profile['n_x'])]
        self.profile['zeta_s'] = profile[(mark + 2*self.profile['n_x']):(mark + 3*self.profile['n_x'])]
        self.profile['kappa_s'] = profile[(mark + 3*self.profile['n_x']):(mark + 4*self.profile['n_x'])]
        self.profile['drmaj_s'] = profile[(mark + 4*self.profile['n_x']):(mark + 5*self.profile['n_x'])]
        self.profile['shat_s'] = profile[(mark + 5*self.profile['n_x']):(mark + 6*self.profile['n_x'])]
        self.profile['s_delta_s'] = profile[(mark + 6*self.profile['n_x']):(mark + 7*self.profile['n_x'])]
        self.profile['s_zeta_z'] = profile[(mark + 7*self.profile['n_x']):(mark + 8*self.profile['n_x'])]
        self.profile['s_kappa_s'] = profile[(mark + 8*self.profile['n_x']):(mark + 9*self.profile['n_x'])]
        self.profile['zmag_s'] = profile[(mark + 9*self.profile['n_x']):(mark + 10*self.profile['n_x'])]
        self.profile['dzmag_s'] = profile[(mark + 10*self.profile['n_x']):(mark + 11*self.profile['n_x'])]
        self.profile['beta_unit_s'] = profile[(mark + 11*self.profile['n_x']):(mark + 12*self.profile['n_x'])]
        self.profile['gamma_e_s'] = profile[(mark + 12*self.profile['n_x']):(mark + 13*self.profile['n_x'])]
        self.profile['gamma_p_s'] = profile[(mark + 13*self.profile['n_x']):(mark + 14*self.profile['n_x'])]
        self.profile['mach_s'] = profile[(mark + 14*self.profile['n_x']):(mark + 15*self.profile['n_x'])]
        self.profile['b_unit_s'] = profile[(mark + 15*self.profile['n_x']):(mark + 16*self.profile['n_x'])]
        self.profile['dr_eodr'] = profile[(mark + 16*self.profile['n_x']):(mark + 17*self.profile['n_x'])]
        self.profile['z_eff_s'] = profile[(mark + 17*self.profile['n_x']):(mark + 18*self.profile['n_x'])]
        self.profile['nu_s'] = profile[(mark + 18*self.profile['n_x']):(mark + 19*self.profile['n_x'])]
        self.profile['w0_s'] = profile[(mark + 19*self.profile['n_x']):(mark + 20*self.profile['n_x'])]
        mark = mark + 20*self.profile['n_x']
        self.profile['box_multiplier'] = profile[mark]
        self.profile['lambda'] = profile[(mark + 1):(mark + 1 + self.profile['n_lambda'])]
        mark = mark + 1 + self.profile['n_lambda']
        self.profile['energy'] = profile[mark:(mark + self.profile['n_energy'])]
        self.profile['lambda_tp'] = profile[mark + self.profile['n_energy']]
        mark = mark + self.profile['n_energy'] + 1
        self.profile['kt_rho'] = profile[mark:(mark + self.profile['n_n'])]
        self.profile['rho_s'] = profile[mark + self.profile['n_n']]
        mark = mark + self.profile['n_n'] + 1
        self.profile['z'] = profile[mark:(mark + self.profile['n_spec'])]
        self.profile['n_fine'] = profile[mark + self.profile['n_spec']]
        self.profile['n_moment'] = profile[mark + self.profile['n_spec'] + 1]
        if int(self.profile['n_theta_plot']) == 1:
            self.profile['theta_plot'] = 0
        else:
            temp = []
            for k in range(int(self.profile['n_theta_plot'])):
                temp.append(-math.pi + 2*math.pi*k/float(self.profile['n_theta_plot']))
            self.profile['theta_plot'] = temp
                

    def read_geometry(self):
        """Reads in geometry_array data.  Output is dictionary of numpy arrays
        with dimensions: n_fine x n_x."""

        import numpy as np
        
        geometry = self.read_file('geometry_arrays', 16)
        if len(geometry) > 0:
            temp = geometry.reshape( (11, self.profile['n_fine'], self.profile['n_x']), order='F')
            self.geometry['v'] = temp[0, :, :]
            self.geometry['gsin'] = temp[1, :, :]
            self.geometry['gcos1'] = temp[2, :, :]
            self.geometry['gcos2'] = temp[3, :, :]
            self.geometry['usin'] = temp[4, :, :]
            self.geometry['ucos'] = temp[5, :, :]
            self.geometry['B'] = temp[6, :, :]
            self.geometry['G_theta'] = temp[7, :, :]
            self.geometry['grad_r'] = temp[8, :, :]
            self.geometry['G_q'] = temp[9, :, :]
            self.geometry['THETA'] = temp[10, :, :]

    def read_freq(self):
        """Reads in frequency data.  Output is dictionary of numpy arrays with
        dimensions: n_n x n_time"""

        import numpy as np

        freq = np.loadtxt(file(self.directory_name + '/freq.out'))
        temp = freq.reshape( (4, self.profile['n_n'], self.t['n_time']), order='F')
        self.freq['(a/c_x)w_{R,n}'] = temp[0, :, :]
        self.freq['(a/c_x)gamma_n'] = temp[1, :, :]
        self.freq['errorin(a/c_x)w_{R,n}'] = temp[2, :, :]
        self.freq['errorin(a/c_x)gamma_n'] = temp[3, :, :]
        self.loaded.append('freq')

    def read_gbflux_i(self):
        """Reads in gbflux_i data.  Output is numpy array with dimensions:
        n_kinetic x n_field x 4 x n_x x n_time"""

        import numpy as np

        gbflux_i = self.read_file('gbflux_i', 12)
        if len(gbflux_i) > 0:
            t = len(gbflux_i)/(self.profile['n_kinetic']*self.profile['n_field']*4*self.profile['n_x'])
            self.gbflux_i = gbflux_i.reshape( (self.profile['n_kinetic'], self.profile['n_field'], 4, self.profile['n_x'], t), order='F')
            self.loaded.append('gbflux_i')

    def read_gbflux_n(self):
        """Reads gbflux_n data.  Output is numpy array with dimensions:
        n_kinetic x n_field x 4 x n_n x n_time"""

        import numpy as np

        gbflux_n = self.read_file('gbflux_n', 12)
        if len(gbflux_n) > 0:
            t = len(gbflux_n)/(self.profile['n_kinetic']*self.profile['n_field']*4*self.profile['n_n'])
            self.gbflux_n = gbflux_n.reshape( (self.profile['n_kinetic'], self.profile['n_field'], 4, self.profile['n_n'], t), order='F')
            self.loaded.append('gbflux_n')

    def read_moment_u(self):
        """Reads in moment_u data.  Output is numpy array with dimensions:
        2 x n_theta_plot x n_x x n_field x n_n x n_time"""

        import numpy as np
        import time
        import fileupload

        starttime = time.time()
        moment_u = self.read_file('moment_u', 12)
        if len(moment_u) > 0:
            t = len(moment_u)/(2*self.profile['n_theta_plot']*self.profile['n_x']*self.profile['n_field']*self.profile['n_n'])
            self.moment_u = moment_u.reshape( (2, self.profile['n_theta_plot'], self.profile['n_x'], self.profile['n_field'], self.profile['n_n'], t), order='F')
            endtime = time.time()
            print "Time: " + str(endtime - starttime)
            self.loaded.append('moment_u')

    def read_moment_n(self):
        """Reads in moment_n data.  Output is numpy array with dimensions:
        2 x n_theta_plot x n_x x n_kinetic x n_n x n_time"""

        import numpy as np
        import time
        import fileupload

        starttime = time.time()
        moment_n = self.read_file('moment_n', 12)
        if len(moment_n) > 0:
            t = len(moment_n)/(2*self.profile['n_theta_plot']*self.profile['n_x']*self.profile['n_kinetic']*self.profile['n_n'])
            self.moment_n = moment_n.reshape( (2, self.profile['n_theta_plot'], self.profile['n_x'], self.profile['n_kinetic'], self.profile['n_n'], t), order='F')
            endtime = time.time()
            print "Time: " + str(endtime - starttime)
            self.loaded.append('moment_n')

    def read_moment_e(self):
        """Reads in moment_e data.  Output is numpy array with dimensions:
        2 x n_theta_plot x n_x x n_kinetic x n_n x n_time"""

        import numpy as np
        import time
        import fileupload

        starttime = time.time()
        moment_e = self.read_file('moment_e', 12)
        if len(moment_e) > 0:
            t = len(moment_e)/(2*self.profile['n_theta_plot']*self.profile['n_x']*self.profile['n_kinetic']*self.profile['n_n'])
            self.moment_e = moment_e.reshape( (2, self.profile['n_theta_plot'], self.profile['n_x'], self.profile['n_kinetic'], self.profile['n_n'], t), order='F')
            endtime = time.time()
            print "Time: " + str(endtime - starttime)
            self.loaded.append('moment_e')

    def read_moment_v(self):
        """Reads in moment_v data.  Output is numpy array with dimensions:
        2 x n_theta_plot x n_x x n_kinetic x n_n x n_time"""

        import numpy as np
        import time
        import fileupload

        starttime = time.time()
        moment_v = self.read_file('moment_v', 12)
        if len(moment_v) > 0:
            midtime = time.time()
            t = len(moment_v)/(2*self.profile['n_theta_plot']*self.profile['n_x']*self.profile['n_kinetic']*self.profile['n_n'])
            self.moment_v = moment_v.reshape( (2, self.profile['n_theta_plot'], self.profile['n_x'], self.profile['n_kinetic'], self.profile['n_n'], t), order='F')
            endtime = time.time()
            print "Time to load data: " + str(midtime - starttime)
            print "Time to reshape data: " + str(endtime - midtime)
            print "Total time: " + str(endtime - starttime)
            self.loaded.append('moment_v')

    def read_moment_zero(self):
        """Reads in moment_zero data.  Output is numpy array with dimensions:
        n_x x n_kinetic x n_moment x n_time"""

        import numpy as np
        import time

        starttime = time.time()
        moment_zero = self.read_file('moment_zero', 12)
        if len(moment_zero) > 0:
            t = len(moment_zero)/(self.profile['n_x']*self.profile['n_kinetic']*self.profile['n_moment'])
            self.moment_zero = moment_zero.reshape( (self.profile['n_x'], self.profile['n_kinetic'], self.profile['n_moment'], t), order='F')
            endtime = time.time()
            print "Time: " + str(endtime - starttime)
            self.loaded.append('moment_zero')

    def read_flux_velocity(self):
        """Reads in flux_velocity data.  Output is numpy array with dimensions:
        n_energy x n_lambda x n_kinetic x n_field x 2 x n_n x n_time"""

        import numpy as np

        flux_velocity = self.read_file('flux_velocity', 12)
        if len(flux_velocity) > 0:
            t = len(flux_velocity)/(self.profile['n_energy']*self.profile['n_lambda']*self.profile['n_kinetic']*self.profile['n_field']*2*self.profile['n_n'])
            self.flux_velocity = flux_velocity.reshape( (self.profile['n_energy'], self.profile['n_lambda'], self.profile['n_kinetic'], self.profile['n_field'], 2, self.profile['n_n'], t), order='F')
            self.loaded.append('flux_velocity')

    def read_k_perp_squared(self):
        """Reads in k_perp_squared data.  Output is numpy array with dimensions:
        n_n x n_time"""

        import numpy as np

        k_perp_squared = self.read_file('k_perp_squared', 12)
        if len(k_perp_squared) > 0:
            t = len(k_perp_squared)/self.profile['n_n']
            self.k_perp_squared = k_perp_squared.reshape( (self.profile['n_n'], t), order='F')
            self.loaded.append(self.k_perp_squared)

    #-----------------------------------#
    # Create data from other previously imported data

    def make_gbflux(self):
        """Makes gbflux.  Output is numpy array with dimensions:
        n_kinetic x n_field x 4 x n_time"""

        import numpy as np

        self.gbflux = np.mean(self.gbflux_i, axis=3)

    def make_diff(self):
        """Makes diff.  Output is dictionary of numpy arrays with
        dimensions: n_kinetic x n_field x n_time"""

        import numpy as np

        if self.gbflux == []:
            self.make_gbflux()
        temp1 = []
        temp2 = []
        temp = []
        for i in range(int(self.profile['n_kinetic'])):
            temp1.append(self.gbflux[i, :, 0, :]/(np.mean(self.profile['dlnndr_s'], axis=1)[i]*np.mean(self.profile['den_s'], axis=1)[i]))
            temp2.append(self.gbflux[i, :, 1, :]/(np.mean(self.profile['dlntdr_s'], axis=1)[i]*np.mean(self.profile['tem_s'], axis=1)[i]*np.mean(self.profile['den_s'], axis=1)[i]))
        temp.append(temp1)
        temp.append(temp2)
        temp = np.swapaxes(np.swapaxes(temp, 0, 1), 1, 2)
        self.diff = np.array(temp)

    def make_diff_i(self):
        """Makes diff_i.  Output is dictionary of numpy arrays with
        dimensions: n_kinetic x n_field x n_x x n_time"""

        import numpy as np

        if self.gbflux_i == []:
            self.read_gbflux_i()
        temp1 = []
        temp2 = []
        temp = []
        for i in range(int(self.profile['n_kinetic'])):
            temp1.append(self.gbflux_i[i, :, 0, :, :]/(np.mean(self.profile['dlnndr_s'], axis=1)[i]*np.mean(self.profile['den_s'], axis=1)[i]))
            temp2.append(self.gbflux_i[i, :, 0, :, :]/(np.mean(self.profile['dlntdr_s'], axis=1)[i]*np.mean(self.profile['tem_s'], axis=1)[i]))
        temp.append(temp1)
        temp.append(temp2)
        temp = np.swapaxes(np.swapaxes(temp, 0, 1), 1, 2)
        self.diff_i = np.array(temp)

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
        """Equilizes the lengths of the different data arrays, in the case that
        the time axis is longer in some than in others."""

        import numpy as np

        temp = []
        for item in self.loaded:
            if item == 't':
                temp.append(self.t['n_time'])
            elif item == 'freq':
                temp.append(len(self.freq.T))
            elif item == 'gbflux_i':
                temp.append(len(self.gbflux_i.T))
            elif item == 'gbflux_n':
                temp.append(len(self.gbflux_n.T))
            elif item == 'moment_u':
                temp.append(len(self.moment_u.T))
            elif item == 'moment_n':
                temp.append(len(self.moment_n.T))
            elif item == 'moment_e':
                temp.append(len(self.moment_e.T))
            elif item == 'moment_v':
                temp.append(len(self.moment_v.T))
            elif item == 'moment_zero':
                temp.append(len(self.moment_zero.T))
            elif item == 'flux_velocity':
                temp.append(len(self.flux_velocity.T))
            elif item == 'k_perp_squared':
                temp.append(len(self.k_perp_squared.T))
        cutoff = min(temp)
        for item in self.loaded:
            if item == 't':
                self.t['t/deltat'] = np.delete(self.t['t/deltat'],
                                               np.s_[cutoff:self.t['n_time']:1],
                                               axis=-1)
                self.t['(cbar_s/a)t'] = np.delete(self.t['(cbar_s/a)t'],
                                                  self.t['n_time']-cutoff,
                                                  axis=-1)
                self.t['n_time'] = cutoff
            elif item == 'freq':
                self.freq['(a/c_x)w_{R,n}'] = np.delete(self.freq['(a/c_x)w_{R,n}'], np.s_[cutoff:len(self.freq['(a/c_x)w_{R,n}'].T):1], axis=-1)
                self.freq['(a/c_x)gamma_n'] = temp[1, :, :]
                self.freq['errorin(a/c_x)w_{R,n}'] = temp[2, :, :]
                self.freq['errorin(a/c_x)gamma_n'] = temp[3, :, :]
            elif item == 'gbflux_i':
                self.gbflux_i = np.delete(self.gbflux_i,
                                          np.s_[cutoff:len(self.gbflux_i.T):1],
                                          axis=-1)
            elif item == 'gbflux_n':
                self.gbflux_n = np.delete(self.gbflux_n,
                                          np.s_[cutoff:len(self.gbflux_n.T):1],
                                          axis=-1)
            elif item == 'moment_u':
                self.moment_u = np.delete(self.moment_u,
                                          np.s_[cutoff:len(self.moment_u.T):1],
                                          axis=-1)
            elif item == 'moment_n':
                self.moment_n = np.delete(self.moment_n,
                                          np.s_[cutoff:len(self.moment_n.T):1],
                                          axis=-1)
            elif item == 'moment_e':
                self.moment_e = np.delete(self.moment_e,
                                          np.s_[cutoff:len(self.moment_e.T):1],
                                                axis=-1)
            elif item == 'moment_v':
                self.moment_v = np.delete(self.moment_v,
                                          np.s_[cutoff:len(self.moment_v.T):1],
                                                axis=-1)
            elif item == 'moment_zero':
                self.moment_zero = np.delete(self.moment_zero,
                                       np.s_[cutoff:len(self.moment_zero.T):1],
                                             axis=-1)
            elif item == 'flux_velocity':
                self.flux_velocity = np.delete(self.flux_velocity,
                                     np.s_[cutoff:len(self.flux_velocity.T):1],
                                               axis=-1)
            elif item == 'k_perp_squared':
                self.k_perp_squared = np.delete(self.k_perp_squared,
                                    np.s_[cutoff:len(self.k_perp_squared.T):1],
                                                axis=-1)
