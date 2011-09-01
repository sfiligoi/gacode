class TGYROData:
    """A class of TGYRO output data.

     Data:

     loc_n_ion
     tgyro_mode
     n_iterations
     n_fields
     n_radial
     dirname = ""
     data

     Example Usage:
         >>> from matplotlib import pyplot
         >>> from pyrats.tgyro.data import TGYROData
         >>> sim1 = TGYROData('$GACODE_ROOT/tgyro/tools/input/treg01')
         >>> pyplot.plot(sim1.get_r(), sim1.get_Te())
         >>> pyplot.show()

"""

    # Methods
    def __init__(self, sim_directory):
        """Constructor reads in data from sim_directory and creates new object.
"""
        self.set_directory(sim_directory)
        self.init_data()
        self.read_data()

    def init_data(self):
        """Initialize object data."""

        self.loc_n_ion    = 0
        self.tgyro_mode   = 0
        self.n_iterations = 0
        self.n_fields     = 0
        self.n_radial     = 0
        self.data         = {}

    def set_directory(self, sim_directory):
        """Set the simulation directory."""

        from os.path import expanduser, expandvars
        path = sim_directory
        self.dirname = expanduser(expandvars(path))

    def get_input(self, input_name):
        """Return the specified variable from input.tgyro.gen.

        input_name  -  requested input

        Ex:    get_input("TGYRO_MODE")
        """

        input_file = file(self.dirname + '/input.tgyro.gen', 'r')
        for line in input_file:
            try:
                if line.split()[1] == input_name:
                    return float(line.split()[0])
            except IndexError:
                print "Cannot find specified input parameter: ", input_name
                return 0

    def read_tgyro_mode(self):
        """Read TGYRO_MODE and store as self.tgyro_mode."""
        
        self.tgyro_mode = self.get_input("TGYRO_MODE")

    def read_num_ions(self):
        """Read LOC_N_ION and store as self.loc_n_ion."""

        self.loc_n_ion = self.get_input("LOC_N_ION")

    def read_data(self):
        """Read in object data."""
        self.read_tgyro_mode()
        if self.tgyro_mode == 2:
            self.read_stabilities()
        else:
            self.read_num_ions()
            self.read_control()
            self.read_chi_e()
            self.read_chi_i(self.loc_n_ion)
            self.read_gyrobohm()
            self.read_profile()
            self.read_geometry()
            self.read_flux(self.loc_n_ion)
            self.read_mflux(self.loc_n_ion)
            self.read_gradient()
            self.read_residual()
            if self.loc_n_ion > 1:
                self.read_profile2()
            if self.loc_n_ion > 2:
                self.read_profile3()
            if self.loc_n_ion > 3:
                self.read_profile4()
            if self.loc_n_ion > 4:
                self.read_profile5()

    def read_control(self):
        """Read control.out to set resolutions."""
        from numpy import loadtxt
        control_file = self.dirname + '/control.out'
        control = loadtxt(file(control_file))
        self.n_radial = int(control[0])
        self.n_fields = int(control[1])
        self.n_iterations = int(control[2])

    def read_file(self, file_name):
        """Read TGYRO output file. Output is data['column_header'][iteration]"""
        from numpy import array

        current_line_number = 0
        elements = {}
        temp = []
        fname = self.dirname + '/' + file_name + '.out'
        raw_data = open(fname, 'r').readlines()
        for line in raw_data:
            if len(line.strip()) > 0:
                if line.strip()[0].isdigit():
                    temp.append(map(float, line.split()))
        data = array(temp)

        # This catches error from columns running into each other
        keywords = raw_data[0].replace('target', 'target ').split()
        for key in keywords:
            elements[key] = range(self.n_iterations + 1)
        for iteration in range(self.n_iterations + 1):
            column = 0
            for key in keywords:
                elements[key][iteration] = data[0:self.n_radial, column]
                column = column + 1
            data = data[self.n_radial:, :]
        return elements

    def read_residual(self):
        """Read residual.out"""
        from numpy import zeros, array
        residual_file=open(self.dirname + '/residual.out', 'r')
        lines = residual_file.readlines()
        residual_file.close()
        count = 0
        iteration = -1
        data = []
        num_fields = len(lines[-1].split())
        local_res = []
        global_res = []
        flux_count = []

        for line in lines:
            line = line.replace(']',' ').replace('[',' ')
            if count % self.n_radial == 0:
                global_res.append(eval(line.split()[3]))
                flux_count.append(eval(line.split()[4]))
                iteration = iteration + 1
                if data:
                    local_res.append(array(data))
                data = [zeros(num_fields)]

            else:
                b = []
                for i in line.split():
                    b.append(eval(i))
                data.append(b)
            count = count + 1
        local_res.append(array(data))
        self.data['local_res'] = local_res
        self.data['global_res'] = global_res
        self.data['flux_count'] = flux_count

    def read_stab_file(self, file_name):
        """Read files generated with stability analysis mode.

        file_name - stability file, eg "wi_elec.out" (Default)

        Call with read_stab_file("wr_elec.out").

        Returns [r, ks, freq]

        """

        from numpy import loadtxt

        full_path = self.dirname + '/' + file_name
        raw_data = loadtxt(file(full_path), skiprows=2, unpack=True)
        r = raw_data[0,:]
        freqs = raw_data[1:,:]
        num_ks = freqs.shape[1]
        ks = loadtxt(file(full_path), skiprows=1, usecols=range(1,num_ks+1))[0]

        return [r, ks, freqs]

    def read_stabilities(self):
        """Read output files from TGYRO_METHOD=2, store into variables.

        Files read:
            wr_ion.out   ->    self.wr_ion
            wi_ion.out   ->    self.wi_ion
            wr_elec.out  ->    self.wr_elec
            wi_elec.out  ->    self.wi_elec

        """

        self.data['wr_ion'] = self.read_stab_file("wr_ion.out")
        self.data['wi_ion'] = self.read_stab_file("wi_ion.out")
        self.data['wr_elec'] = self.read_stab_file("wr_elec.out")
        self.data['wi_elec'] = self.read_stab_file("wi_elec.out")

        self.data['r/a'] = self.data['wr_ion'][0]
        self.n_radial = len(self.data['r/a'])

    def read_flux(self, num_ions = 1):
        """Read flux_e.out, flux_i(2-5).out, flux_target.out."""
        self.data.update(self.read_file('flux_e'))
        self.data.update(self.read_file('flux_i'))
        if num_ions > 1:
           self.data.update(self.read_file('flux_i2'))
        if num_ions > 2:
           self.data.update(self.read_file('flux_i3'))
        if num_ions > 3:
           self.data.update(self.read_file('flux_i4'))
        if num_ions > 4:
           self.data.update(self.read_file('flux_i5'))
        if num_ions > 5:
            print "Too many ions: ", num_ions
            print "Only the first 5 will be read."
        self.data.update(self.read_file('flux_target'))

    def read_mflux(self, num_ions = 1):
        """Read mflux_e.out, mflux_i(2-5).out, mflux_target.out."""
        self.data.update(self.read_file('mflux_e'))
        self.data.update(self.read_file('mflux_i'))
        if num_ions > 1:
           self.data.update(self.read_file('mflux_i2'))
        if num_ions > 2:
           self.data.update(self.read_file('mflux_i3'))
        if num_ions > 3:
           self.data.update(self.read_file('mflux_i4'))
        if num_ions > 4:
           self.data.update(self.read_file('mflux_i5'))
        if num_ions > 5:
            print "Too many ions: ", num_ions
            print "Only the first 5 will be read."
        self.data.update(self.read_file('mflux_target'))

    def read_chi_e(self):
        """Read in chi_e.out and store in self.chi_e."""
        self.data.update(self.read_file('chi_e'))

    def read_chi_i(self, num_ions = 1):
        """Read in chi_i.out and store in self.chi_i."""
        self.data.update(self.read_file('chi_i'))
        if num_ions > 1:
            self.data.update(self.read_file('chi_i2'))
        if num_ions > 2:
            self.data.update(self.read_file('chi_i3'))
        if num_ions > 3:
            self.data.update(self.read_file('chi_i4'))
        if num_ions > 4:
            self.data.update(self.read_file('chi_i5'))
        if num_ions > 5:
            print "Strange number of ions: ", num_ions
            print "Only the first 5 will be read."

    def read_gyrobohm(self):
        """Read and store gyrobohm.out in self.gyro_bohm_unit."""
        self.data.update(self.read_file('gyrobohm'))

    def read_profile(self):
        """Read and store profile.out in self.profile."""
        self.data.update(self.read_file('profile'))

    def read_profile2(self):
        """Read and store profile2.out in self.profile2."""
        self.data.update(self.read_file('profile2'))

    def read_profile3(self):
        """Read and store profile3.out in self.profile3."""
        self.data.update(self.read_file('profile3'))

    def read_profile4(self):
        """Read and store profile4.out in self.profile4."""
        self.data.update(self.read_file('profile4'))

    def read_profile5(self):
        """Read and store profile5.out in self.profile5."""
        self.data.update(self.read_file('profile5'))

    def read_geometry(self):
        """Read and store geometry.out in self.geometry."""
        self.data.update(self.read_file('geometry'))

    def read_gradient(self):
        """Read and store gradient.out in self.gradient."""
        self.data.update(self.read_file('gradient'))


    
    # ----------------------------------------- #
    # Get data back

    def get_local_res(self, field=1, iteration=-1):
        try:
            return self.data['local_res'][iteration][:,2*field-1]
        except IndexError:
            print "Invalid field residual."
            print self.data['local_res']
            return None

    def get_stability_at_radius(self, radius=0, frequency='r', direction='ion'):
        """Get frequency vs. ky at specified radius from stability analysis.

        Parameters:

            radius        -       index of requested radius (Default: 0)
                                  integer ranging from 0 : n_r-1

            frequency     -       real or imaginary spectrum
                                  'r'     : real (Default)
                                  'i'     : imaginary

            direction     -       direction of spectrum
                                  'ion'   : ion direction (Default)
                                  'elec'  : electron direction

        Usage:

            [ky, wi_ion] = get_stability_at_radius(0, 'i', 'ion')

        """

        # Check for valid inputs


        if radius < 0 or radius > self.n_radial-1:
            print "Error in get_stability_at_radius. Invalid radial point: ", radius
            print "Valid range: 0 - ", self.n_radial-1
            return

        if frequency is not 'r' and frequency is not 'i':
            print "Error in get_stability_at_radius. Invalid frequency: ", frequency
            print "Valid options are 'r' and 'i'"
            return

        if direction is not 'ion' and direction is not 'elec':
            print "Error in get_stability_at_radius. Invalid frequency: ", frequency
            print "Valid options are 'ion' and 'elec'"
            return

        # Return requested data
        if direction is 'ion':
            if frequency is 'r':
                k = self.data['wr_ion'][1]
                f = self.data['wr_ion'][2][:, radius]
            else:
                k = self.data['wi_ion'][1]
                f = self.data['wi_ion'][2][:, radius]
        else:
            if frequency is 'r':
                k = self.data['wr_elec'][1]
                f = self.data['wr_elec'][2][:, radius]
            else:
                k = self.data['wi_elec'][1]
                f = self.data['wi_elec'][2][:, radius]

        return [k, f]

    def get_most_unstable_at_radius(self, radius=0, direction='ion'):
        """Return the most unstable mode at specified radius in the specified
        direction.

        Parameters:

           radius     -   requested radial index (Default: 0)
                          integer between 0 and n_r-1

           direction  -   'ion'  most unstable ion mode
                          'elec' most unstable electron mode

        Usage:

           [k_e, omega_e, gamma_e] = get_most_unstable_at_radius(0, 'elec')

        """

        # Check input sanity

        if radius < 0 or radius > self.n_radial-1:
            print "Error in get_most_unstable_at_radius. Invalid radial point: ", radius
            print "Valid range: 0 - ", self.n_radial-1
            return

        if direction is not 'ion' and direction is not 'elec':
            print "Error in get_most_unstable_at_radius. Invalid frequency: ", frequency
            print "Valid options are 'ion' and 'elec'"
            return

        # Get ion and electron data at requested radius

        ks, gammas = self.get_stability_at_radius(radius, 'i', direction)
        k2, omegas = self.get_stability_at_radius(radius, 'r', direction)

        gamma_max = max(gammas)
        max_index = gammas.argmax()
        k_max = ks[max_index]
        omega_max = omegas[max_index]
        
        return k_max, omega_max, gamma_max

# Analysis Tools

    def make_gradient_vs_field_space(self, r=1, evolve_field=1, grad='a/LTi', profile='ti'):
        """Return matrix of Space[iteration][gradient,profile,residual] at given radial point.

        eg: Space[[4.0, 1.7, 2.5], [3.6, 1.5, 0.1]] if grad 4.0->3.6 while Ti->1.5
            and local residual went from 2.5 -> 0.1

        """
        import numpy
        space=[]
        for iteration in range(self.n_iterations + 1):
            #print "iteration", iteration, "r", r
            loc_grad = self.data[grad][iteration][r]
            loc_prof = self.data[profile][iteration][r]
            loc_resi = self.get_local_res(evolve_field, iteration)[r]
            #print loc_grad, loc_prof, loc_resi
            space.append([loc_grad, loc_prof, loc_resi])

        return numpy.array(space)             

    def make_gradient_vs_field_plot(self, r=1, evolve_field=1, grad='a/LTi', profile='ti', arrows=False):
        """Return a pyplot figure of gradient_vs_field_space."""
        import matplotlib.pyplot as plt
        from pylab import Arrow
        space = self.make_gradient_vs_field_space(r, evolve_field, grad, profile)
        fig = plt.figure()
        walk = fig.add_subplot(111)
        x=space[:,0]
        y=space[:,1]
        z=space[:,2]
        cax = walk.scatter(x,y,c=z,s=40,vmin=min(z),vmax=max(z))
        walk.set_xlabel(grad)
        walk.set_ylabel(profile)
        cbar = plt.colorbar(cax)
        cbar.set_label("Local Residual")
        walk.set_title("Radial Point " + str(r))
        xmin, xmax = plt.xlim()
        ymin, ymax = plt.ylim()
        deltax = xmax - xmin
        deltay = ymax - ymin

        if arrows:
            arrow_width = 0.05 * deltax
            for point in range(len(x)-1):
                dx = x[point+1] - x[point]
                dy = y[point+1] - y[point]
                arrow = Arrow(x[point], y[point], dx, dy, alpha=0.75, width=arrow_width, \
                              linewidth=0, edgecolor='none', facecolor='black')
                walk.add_patch(arrow)

        for i in range(len(x)-1):
            if i % 5 == 0:
                walk.annotate(str(i), (x[i], y[i]), xycoords='data', size=15)           
        walk.annotate("Start", (x[0], y[0]), xycoords='data', xytext=((xmax - 0.3 * deltax), (ymin + 0.05 * deltay)), \
                      textcoords='data', size=15, \
                      arrowprops=dict(width=3, headwidth=5, shrink=0.05),)
        walk.annotate("End, Iteration "+ str(len(x)-1), (x[-1], y[-1]), xycoords='data', xytext=((xmin + 0.05 * deltax), \
                      (ymax - 0.1 * deltay)), \
                      textcoords='data', size=15, \
                      arrowprops=dict(width=3, headwidth=5, shrink=0.05, facecolor="red"),)
        minpos = 0
        for j in range(len(z-1)):
            if z[j] == min(z):
                minpos = j
        walk.annotate("Min Residual = " + str(z[minpos]), (x[minpos], y[minpos]), xycoords='data', \
                      xytext=((xmin + 0.05 * deltax), (ymin + 0.05 * deltay)), \
                      textcoords='data', size=15, arrowprops=dict(width=3, headwidth=5, shrink=0.05, facecolor="black"))

        return fig
