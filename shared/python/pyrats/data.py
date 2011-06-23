"""data.py contains data classes for pyrats.

    Contents:
        TGYROData

"""

class TGYROData:
    """A class of TGYRO output data.

     Data:

     n_iterations
     n_fields
     n_radial
     directory_name = ""
     chi_e = []
     chi_i = []
     gyro_bohm_unit = []
     profile = []
     geometry = []
     flux_e = []
     flux_i = []
     flux_target = []
     mflux_e = []
     mflux_i = []
     mflux_target = []
     gradient = []
     local_res = []
     global_res = []
     wr_ion = []
     wi_ion = []
     wr_elec = []
     wi_elec = []
     r
     (Optional)
     chi_i2 = []
     chi_i3 = []
     chi_i4 = []
     chi_i5 = []
     flux_i2 = []
     flux_i3 = []
     flux_i4 = []
     flux_i5 = []
     mflux_i2 = []
     mflux_i3 = []
     mflux_i4 = []
     mflux_i5 = []

     Example Usage:
         >>>from matplotlib import pyplot
         >>>from pyrats.data import TGYROData
         >>>sim1 = TGYROData('$GACODE_ROOT/tgyro/tools/input/treg01')
         >>>pyplot.plot(sim1.get_r(), sim1.get_Te())
         >>>pyplot.show()

"""

    # Methods
    def __init__(self, sim_directory):
        """Constructor reads in data from sim_directory and creates new object."""
        self.set_directory(sim_directory)
        self.init_data()
        self.read_data()

    def init_data(self):
        """Initialize object data."""

        self.loc_n_ion = 0
        self.read_num_ions()
        self.tgyro_mode = 0
        self.n_iterations = 0
        self.n_fields = 0
        self.n_radial = 0
        self.chi_e = []
        self.chi_i = []
        self.gyro_bohm_unit = []
        self.profile = []
        self.geometry = []
        self.flux_e = []
        self.flux_i = []
        self.flux_target = []
        self.mflux_e = []
        self.mflux_i = []
        self.mflux_target = []
        self.gradient = []
        self.local_res = []
        self.global_res = []
        self.flux_count = []
        self.wr_ion = []
        self.wi_ion = []
        self.wr_elec = []
        self.wi_elec = []
        self.r = 0
        if self.loc_n_ion > 1:
            self.chi_i2 = []
            self.flux_i2 = []
            self.mflux_i2 = []
        if self.loc_n_ion > 2:
            self.chi_i3 = []
            self.flux_i3 = []
            self.mflux_i3 = []
        if self.loc_n_ion > 3:
            self.chi_i4 = []
            self.flux_i4 = []
            self.mflux_i4 = []
        if self.loc_n_ion > 4:
            self.chi_i5 = []
            self.flux_i5 = []
            self.mflux_i5 = []

    def set_directory(self, sim_directory):
        """Set the simulation directory."""

        from os.path import expanduser, expandvars
        path = sim_directory
        self.directory_name = expanduser(expandvars(path))

    def tgyro_get_input(self, input_name):
        """Return the specified variable from input.tgyro.gen.

        input_name  -  requested input

        Ex:    tgyro_get_input("TGYRO_MODE")
        """

        input_file = file(self.directory_name + '/input.tgyro.gen', 'r')
        for line in input_file:
            try:
                if line.split()[1] == input_name:
                    return float(line.split()[0])
            except IndexError:
                print "Cannot find specified input parameter: ", input_name
                return 0

    def read_tgyro_mode(self):
        """Read TGYRO_MODE and store as self.tgyro_mode."""
        
        self.tgyro_mode = self.tgyro_get_input("TGYRO_MODE")

    def read_num_ions(self):
        """Read LOC_N_ION and store as self.loc_n_ion."""

        self.loc_n_ion = self.tgyro_get_input("LOC_N_ION")

    def read_data(self):
        """Read in object data."""
        self.read_tgyro_mode()
        if self.tgyro_mode == 2:
            self.read_stabilities()
        else:
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
        control_file = self.directory_name + '/control.out'
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
        raw_data = open(file_name, 'r').readlines()
        for line in raw_data:
            if len(line.strip()) > 0:
                if line.strip()[0].isdigit():
                    temp.append(line.split())
        data = array(temp)

        # This catches error from columns running into each other
        keywords = raw_data[0].replace('target', 'target ').split()
        for key in keywords:
            elements[key] = range(self.n_iterations + 1)
        for iteration in range(self.n_iterations + 1):
            column = 0
            for key in keywords:
                elements[key][iteration] = data[1:self.n_radial, column]
                column = column + 1
            data = data[self.n_radial:, :]
        return elements

    def read_residual(self):
        """Read residual.out"""
        from numpy import zeros, array
        residual_file=open(self.directory_name + '/residual.out', 'r')
        lines = residual_file.readlines()
        residual_file.close()
        count = 0
        iteration = -1
        data = []
        num_fields = len(lines[-1].split())

        for line in lines:
            line = line.replace(']',' ').replace('[',' ')
            if count % self.n_radial == 0:
                self.global_res.append(eval(line.split()[3]))
                self.flux_count.append(eval(line.split()[4]))
                iteration = iteration + 1
                if data:
                    self.local_res.append(array(data))
                data = [zeros(num_fields)]

            else:
                b = []
                for i in line.split():
                    b.append(eval(i))
                data.append(b)
            count = count + 1
        self.local_res.append(array(data))    

    def read_stab_file(self, file_name):
        """Read files generated with stability analysis mode.

        file_name - stability file, eg "wi_elec.out" (Default)

        Call with read_stab_file("wr_elec.out").

        Returns [r, ks, freq]

        """

        from numpy import loadtxt

        full_path = self.directory_name + '/' + file_name
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

        self.wr_ion = self.read_stab_file("wr_ion.out")
        self.wi_ion = self.read_stab_file("wi_ion.out")
        self.wr_elec = self.read_stab_file("wr_elec.out")
        self.wi_elec = self.read_stab_file("wi_elec.out")

        self.r = self.wr_ion[0]
        self.n_radial = len(self.r)

    def read_flux(self, num_ions = 1):
        """Read flux_e.out, flux_i(2-5).out, flux_target.out."""
        self.flux_e = self.read_file(self.directory_name + '/flux_e.out')
        self.flux_i = self.read_file(self.directory_name + '/flux_i.out')
        if num_ions > 1:
           self.flux_i2 = self.read_file(self.directory_name + '/flux_i2.out')
        if num_ions > 2:
           self.flux_i3 = self.read_file(self.directory_name + '/flux_i3.out')
        if num_ions > 3:
           self.flux_i4 = self.read_file(self.directory_name + '/flux_i4.out')
        if num_ions > 4:
           self.flux_i5 = self.read_file(self.directory_name + '/flux_i5.out')
        if num_ions > 5:
            print "Too many ions: ", num_ions
            print "Only the first 5 will be read."
        self.flux_target = self.read_file(self.directory_name + '/flux_target.out')

    def read_mflux(self, num_ions = 1):
        """Read mflux_e.out, mflux_i(2-5).out, mflux_target.out."""
        self.mflux_e = self.read_file(self.directory_name + '/mflux_e.out')
        self.mflux_i = self.read_file(self.directory_name + '/mflux_i.out')
        if num_ions > 1:
           self.mflux_i2 = self.read_file(self.directory_name + '/mflux_i2.out')
        if num_ions > 2:
           self.mflux_i3 = self.read_file(self.directory_name + '/mflux_i3.out')
        if num_ions > 3:
           self.mflux_i4 = self.read_file(self.directory_name + '/mflux_i4.out')
        if num_ions > 4:
           self.mflux_i5 = self.read_file(self.directory_name + '/mflux_i5.out')
        if num_ions > 5:
            print "Too many ions: ", num_ions
            print "Only the first 5 will be read."
        self.flux_target = self.read_file(self.directory_name + '/flux_target.out')

    def read_chi_e(self):
        """Read in chi_e.out and store in self.chi_e."""
        self.chi_e = self.read_file(self.directory_name + '/chi_e.out')

    def read_chi_i(self, num_ions = 1):
        """Read in chi_i.out and store in self.chi_i."""
        self.chi_i = self.read_file(self.directory_name + '/chi_i.out')
        if num_ions > 1:
            self.chi_i2 = self.read_file(self.directory_name + '/chi_i2.out')
        if num_ions > 2:
            self.chi_i3 = self.read_file(self.directory_name + '/chi_i3.out')
        if num_ions > 3:
            self.chi_i4 = self.read_file(self.directory_name + '/chi_i4.out')
        if num_ions > 4:
            self.chi_i5 = self.read_file(self.directory_name + '/chi_i5.out')
        if num_ions > 5:
            print "Strange number of ions: ", num_ions
            print "Only the first 5 will be read."

    def read_gyrobohm(self):
        """Read and store gyrobohm.out in self.gyro_bohm_unit."""
        self.gyro_bohm_unit = self.read_file(self.directory_name + '/gyrobohm.out')

    def read_profile(self):
        """Read and store profile.out in self.profile."""
        self.profile = self.read_file(self.directory_name + '/profile.out')

    def read_profile2(self):
        """Read and store profile2.out in self.profile2."""
        self.profile2 = self.read_file(self.directory_name + '/profile2.out')

    def read_profile3(self):
        """Read and store profile3.out in self.profile3."""
        self.profile3 = self.read_file(self.directory_name + '/profile3.out')

    def read_profile4(self):
        """Read and store profile4.out in self.profile4."""
        self.profile4 = self.read_file(self.directory_name + '/profile4.out')

    def read_profile5(self):
        """Read and store profile5.out in self.profile5."""
        self.profile5 = self.read_file(self.directory_name + '/profile.out5')

    def read_geometry(self):
        """Read and store geometry.out in self.geometry."""
        self.geometry = self.read_file(self.directory_name + '/geometry.out')
        self.r = self.geometry['r/a'][-1]

    def read_gradient(self):
        """Read and store gradient.out in self.gradient."""
        self.gradient = self.read_file(self.directory_name + '/gradient.out')


    
    # ----------------------------------------- #
    # Get data back
    def get_r(self):
        """Return r/a."""
        return self.r

    def get_Te(self, iteration=-1):
        """Return Te for specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
        """
        return self.profile['te'][iteration]

    def get_Ti(self, iteration=-1):
        """Return Ti for specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
        """
        return self.profile['ti'][iteration]

    def get_Ti2(self, iteration=-1):
        """Return Ti2 for specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
        """
        return self.profile2['Ti'][iteration]

    def get_Ti3(self, iteration=-1):
        """Return Ti3 for specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
        """
        return self.profile3['Ti'][iteration]

    def get_Ti4(self, iteration=-1):
        """Return Ti4 for specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
        """
        return self.profile4['Ti'][iteration]

    def get_Ti5(self, iteration=-1):
        """Return Ti5 for specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
        """
        return self.profile5['Ti'][iteration]

    def get_ne(self, iteration=-1):
        """Return ne for specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
        """
        return self.profile['ne'][iteration]

    def get_ni(self, iteration=-1):
        """Return ni for specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
        """
        return self.profile['ni'][iteration]

    def get_ni2(self, iteration=-1):
        """Return ni2 for specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
        """
        return self.profile2['ni'][iteration]

    def get_ni3(self, iteration=-1):
        """Return ni3 for specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
        """
        return self.profile3['ni'][iteration]

    def get_ni4(self, iteration=-1):
        """Return ni4 for specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
        """
        return self.profile4['ni'][iteration]

    def get_ni5(self, iteration=-1):
        """Return ni5 for specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
        """
        return self.profile5['ni'][iteration]

    def get_chi_e_turb(self, iteration=-1, mks=False):
        """Return turbulent chi_e from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        chis = self.chi_e['chie_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Chi_GB'][iteration]
            chis = chis * GB_factors

        return chis

    def get_chi_i_turb(self, iteration=-1, mks=False):
        """Return turbulent chi_i from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        chis = self.chi_i['chii_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Chi_GB'][iteration]
            chis = chis * GB_factors

        return chis

    def get_chi_i2_turb(self, iteration=-1, mks=False):
        """Return turbulent chi_i2 from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        chis = self.chi_i2['chii_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Chi_GB'][iteration]
            chis = chis * GB_factors

        return chis

    def get_chi_i3_turb(self, iteration=-1, mks=False):
        """Return turbulent chi_i3 from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        chis = self.chi_i3['chii_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Chi_GB'][iteration]
            chis = chis * GB_factors

        return chis

    def get_chi_i4_turb(self, iteration=-1, mks=False):
        """Return turbulent chi_i4 from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        chis = self.chi_i4['chii_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Chi_GB'][iteration]
            chis = chis * GB_factors

        return chis

    def get_chi_i5_turb(self, iteration=-1, mks=False):
        """Return turbulent chi_i5 from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        chis = self.chi_i5['chii_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Chi_GB'][iteration]
            chis = chis * GB_factors

        return chis

    def get_chi_e_neo(self, iteration=-1, mks=False):
        """Return neoclassical chi_e from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        chis = self.chi_e['chie_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Chi_GB'][iteration]
            chis = chis * GB_factors

        return chis

    def get_chi_i_neo(self, iteration=-1, mks=False):
        """Return neoclassical chi_i from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        chis = self.chi_i['chii_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Chi_GB'][iteration]
            chis = chis * GB_factors

        return chis

    def get_chi_i2_neo(self, iteration=-1, mks=False):
        """Return neoclassical chi_i2 from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        chis = self.chi_i2['chii_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Chi_GB'][iteration]
            chis = chis * GB_factors

        return chis

    def get_chi_i3_neo(self, iteration=-1, mks=False):
        """Return neoclassical chi_i3 from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        chis = self.chi_i3['chii_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Chi_GB'][iteration]
            chis = chis * GB_factors

        return chis

    def get_chi_i4_neo(self, iteration=-1, mks=False):
        """Return neoclassical chi_i4 from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        chis = self.chi_i4['chii_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Chi_GB'][iteration]
            chis = chis * GB_factors

        return chis

    def get_chi_i5_neo(self, iteration=-1, mks=False):
        """Return neoclassical chi_i5 from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        chis = self.chi_i['chii_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Chi_GB'][iteration]
            chis = chis * GB_factors

        return chis

    def get_flux_e_target(self, iteration=-1, mks=False):
        """Return target electron heat flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """

        flux = self.flux_target['eflux_e_target'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_flux_i_target(self, iteration=-1, mks=False):
        """Return target ion heat flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """

        flux = self.flux_target['eflux_i_target'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_flux_e_turb(self, iteration=-1, mks=False):
        """Return turbulent electron heat flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.flux_e['eflux_e_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_flux_i_turb(self, iteration=-1, mks=False):
        """Return turbulent ion heat flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.flux_i['eflux_i_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_flux_i2_turb(self, iteration=-1, mks=False):
        """Return turbulent ion2 heat flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.flux_i2['eflux_i_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_flux_i3_turb(self, iteration=-1, mks=False):
        """Return turbulent ion3 heat flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.flux_i3['eflux_i_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_flux_i4_turb(self, iteration=-1, mks=False):
        """Return turbulent ion4 heat flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.flux_i4['eflux_i_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_flux_i5_turb(self, iteration=-1, mks=False):
        """Return turbulent ion5 heat flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.flux_i5['eflux_i_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_flux_e_neo(self, iteration=-1, mks=False):
        """Return neoclassical electron heat flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.flux_e['eflux_e_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_flux_i_neo(self, iteration=-1, mks=False):
        """Return neoclassical ion heat flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.flux_i['eflux_i_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_flux_i2_neo(self, iteration=-1, mks=False):
        """Return neoclassical ion2 heat flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.flux_i2['eflux_i_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_flux_i3_neo(self, iteration=-1, mks=False):
        """Return neoclassical ion3 heat flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.flux_i3['eflux_i_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_flux_i4_neo(self, iteration=-1, mks=False):
        """Return neoclassical ion4 heat flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.flux_i4['eflux_i_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_flux_i5_neo(self, iteration=-1, mks=False):
        """Return neoclassical ion5 heat flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.flux_i5['eflux_i_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_mflux_target(self, iteration=-1, mks=False):
        """Return target momentum flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """

        flux = self.mflux_target['mflux_target'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_mflux_e_turb(self, iteration=-1, mks=False):
        """Return turbulent electron momentum flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.mflux_e['mflux_e_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_mflux_i_turb(self, iteration=-1, mks=False):
        """Return turbulent ion momentum flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.mflux_i['mflux_i_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_mflux_i2_turb(self, iteration=-1, mks=False):
        """Return turbulent ion2 momentum flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.mflux_i2['mflux_i_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_mflux_i3_turb(self, iteration=-1, mks=False):
        """Return turbulent ion3 momentum flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.mflux_i3['mflux_i_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_mflux_i4_turb(self, iteration=-1, mks=False):
        """Return turbulent ion4 momentum flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.mflux_i4['mflux_i_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_mflux_i5_turb(self, iteration=-1, mks=False):
        """Return turbulent ion5 momentum flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.mflux_i5['mflux_i_tur'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_mflux_e_neo(self, iteration=-1, mks=False):
        """Return neoclassical electron momentum flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.mflux_e['mflux_e_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_mflux_i_neo(self, iteration=-1, mks=False):
        """Return neoclassical ion momentum flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.mflux_i['mflux_i_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_mflux_i2_neo(self, iteration=-1, mks=False):
        """Return neoclassical ion2 momentum flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.mflux_i2['mflux_i_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_mflux_i3_neo(self, iteration=-1, mks=False):
        """Return neoclassical ion3 momentum flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.mflux_i3['mflux_i_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_mflux_i4_neo(self, iteration=-1, mks=False):
        """Return neoclassical ion4 momentum flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.mflux_i4['mflux_i_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_mflux_i5_neo(self, iteration=-1, mks=False):
        """Return neoclassical ion5 momentum flux from specified iteration.
           Keywords:
               iteration           TGYRO iteration to return, default is last
               mks                 whether to return in mks units, default is False
        """
        flux = self.mflux_i5['mflux_i_neo'][iteration]

        if mks:
            GB_factors = self.gyro_bohm_unit['Q_GB'][iteration]
            flux = flux * GB_factors

        return flux

    def get_geometry_factor(self, factor, iteration=-1):
        """Return the specified geometric factor for the given iteration.
           Keywords:
               factor              Geometric Factor, options:
                                       'rho', 'q', 's', 'kappa', 's_kappa', 'delta', 
                                       's_delta', 'shift', 'rmaj/a', 'b_unit'
               iteration           TGYRO iteration to return, default is first
        """
        if factor not in self.geometry.keys():
            print "ERROR: invalid geometric choice. Valid options:"
            print self.geometry.keys()
            return

        return self.geometry[factor][iteration]

    def get_run(self):
        """Return the contents of run.out"""
        run_file = open(self.directory_name + '/run.out', 'r')
        return run_file.read()

    def get_gradient_factor(self, factor, iteration=-1):
        """Return the specified gradient factor for the given iteration.
           Keywords:
               factor              Gradient Factor, options:
                                       'r/a' 'a/Lni' 'a/Lne' 'a/LTi' 'a/LTe' 'a/Lp'
                                       'a*gamma_e/cs' 'a*gamma_p/cs'

               iteration           TGYRO iteration to return, default is last
        """
        if factor not in self.gradient.keys():
            print "ERROR: invalid gradient choice. Valid options:"
            print self.gradient.keys()
            return

        return self.gradient[factor][iteration]

    def get_profile_factor(self, factor, iteration=-1):
        """Return the specified profile factor for the given iteration.
           Keywords:
               factor              Profile Factor, options:
                                       'r/a' 'ni' 'ne' 'ti' 'te' 'ti/te'
                                       'betae_unit' 'M=wR/cs'

               iteration           TGYRO iteration to return, default is last
        """
        if factor not in self.profile.keys():
            print "ERROR: invalid profile choice. Valid options:"
            print self.profile.keys()
            return

        return self.profile[factor][iteration]

    def get_profile2_factor(self, factor, iteration=-1):
        """Return the specified profile2 factor for the given iteration.
           Keywords:
               factor              Profile2 Factor, options:
                                       'r/a' 'ni' 'a/Lni' 'Ti' 'a/LTi'

               iteration           TGYRO iteration to return, default is last
        """
        if factor not in self.profile2.keys():
            print "ERROR: invalid profile2 choice. Valid options:"
            print self.profile2.keys()
            return

        return self.profile2[factor][iteration]

    def get_profile3_factor(self, factor, iteration=-1):
        """Return the specified profile3 factor for the given iteration.
           Keywords:
               factor              Profile3 Factor, options:
                                       'r/a' 'ni' 'a/Lni' 'Ti' 'a/LTi'

               iteration           TGYRO iteration to return, default is last
        """
        if factor not in self.profile3.keys():
            print "ERROR: invalid profile3 choice. Valid options:"
            print self.profile3.keys()
            return

        return self.profile3[factor][iteration]

    def get_profile4_factor(self, factor, iteration=-1):
        """Return the specified profile4 factor for the given iteration.
           Keywords:
               factor              Profile4 Factor, options:
                                       'r/a' 'ni' 'a/Lni' 'Ti' 'a/LTi'

               iteration           TGYRO iteration to return, default is last
        """
        if factor not in self.profile4.keys():
            print "ERROR: invalid profile4 choice. Valid options:"
            print self.profile4.keys()
            return

        return self.profile4[factor][iteration]

    def get_profile5_factor(self, factor, iteration=-1):
        """Return the specified profile5 factor for the given iteration.
           Keywords:
               factor              Profile5 Factor, options:
                                       'r/a' 'ni' 'a/Lni' 'Ti' 'a/LTi'

               iteration           TGYRO iteration to return, default is last
        """
        if factor not in self.profile5.keys():
            print "ERROR: invalid profile5 choice. Valid options:"
            print self.profile5.keys()
            return

        return self.profile5[factor][iteration]

    def get_global_res(self, iteration=-1):
        """Return the global residual from the given iteration."""
        return self.global_res[iteration]

    def get_local_res(self, field=1, iteration=-1):
        try:
            return self.local_res[iteration][:,2*field-1]
        except IndexError:
            print "Invalid field residual."
            print self.local_res
            return None

    def get_flux_count(self, iteration=-1):
        """Get the total number of calls to flux driver up to given iteration."""
        return self.flux_count[iteration]

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
                k = self.wr_ion[1]
                f = self.wr_ion[2][:, radius]
            else:
                k = self.wi_ion[1]
                f = self.wi_ion[2][:, radius]
        else:
            if frequency is 'r':
                k = self.wr_elec[1]
                f = self.wr_elec[2][:, radius]
            else:
                k = self.wi_elec[1]
                f = self.wi_elec[2][:, radius]

        return [k, f]

    def get_most_unstable_at_radius(self, radius=0, direction='ion'):
        """Return the most unstable mode at specified radius in the specified direction.

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
            loc_grad = self.get_gradient_factor(grad, iteration)[r]
            loc_prof = self.get_profile_factor(profile, iteration)[r]
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
