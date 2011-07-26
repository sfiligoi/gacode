"""data.py contains data classes for pyrats.

    Contents:
        TGYROData
        NEOData
        GYROData

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
        """Constructor reads in data from sim_directory and creates new object.
"""
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

    def get_input(self, input_name):
        """Return the specified variable from input.tgyro.gen.

        input_name  -  requested input

        Ex:    get_input("TGYRO_MODE")
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
                elements[key][iteration] = data[0:self.n_radial, column]
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
        self.flux_target = self.read_file(self.directory_name +
                                          '/flux_target.out')

    def read_mflux(self, num_ions = 1):
        """Read mflux_e.out, mflux_i(2-5).out, mflux_target.out."""
        self.mflux_e = self.read_file(self.directory_name + '/mflux_e.out')
        self.mflux_i = self.read_file(self.directory_name + '/mflux_i.out')
        if num_ions > 1:
           self.mflux_i2 = self.read_file(self.directory_name +
                                          '/mflux_i2.out')
        if num_ions > 2:
           self.mflux_i3 = self.read_file(self.directory_name +
                                          '/mflux_i3.out')
        if num_ions > 3:
           self.mflux_i4 = self.read_file(self.directory_name +
                                          '/mflux_i4.out')
        if num_ions > 4:
           self.mflux_i5 = self.read_file(self.directory_name +
                                          '/mflux_i5.out')
        if num_ions > 5:
            print "Too many ions: ", num_ions
            print "Only the first 5 will be read."
        self.flux_target = self.read_file(self.directory_name +
                                          '/flux_target.out')

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
        self.gyro_bohm_unit = self.read_file(self.directory_name +
                                             '/gyrobohm.out')

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
        self.profile5 = self.read_file(self.directory_name + '/profile5.out')

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
        """Get the total number of calls to flux driver up to given iteration.
        """
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

#==============================================================================
#==============================================================================

class NEOOutput:
    """A class which holds individual pieces of NEO output data.

    This class is just a simple holder that can attach the units and a
    description to a piece of data.  It's only purpose is to be used by NEOData.

    Data:
    
    data = {}
    units = ''
    descriptor = ''
"""

    def __init__(self, data, units, descriptor):
        self.data=data
        self.units=units
        self.descriptor=descriptor

class NEOData:
    """A class of NEO output data.

    Data (things that actually get used):

    master = ""
    directory_name = ""
    fignum = 1
    plotcounter = 1
    toplot = []
    transport = {}
    HH_theory = {}
    CH_theory = {}
    TG_theory = {}
    S_theory = {}
    HR_theory = {}
    HS_theory = {}
    control = {}

    Temporary data (data which is stored in one of the above objects):

    n_species = {}
    n_energy = {}
    n_xi = {}
    n_theta = {}
    theta_gpoints = {}
    n_radial = {}
    radial_gpoints = {}
    r = {}
    dphi0dr = {}
    q = {}
    rho_star = {}
    rmaj_over_a = {}
    omega_rot = {}
    omega_rot_deriv = {}
    nnorm = {}
    tnorm = {}
    n0_over_nnorm = {}
    t0_over_tnorm = {}
    a_over_Ln = {}
    a_over_Lt = {}
    inv_tau_self = {}
    f = {}
    SIM_PHI_by_theta = {}
    HH_GAMMA = {}
    HH_Qi = {}
    HH_Qe = {}
    HH_jboot = {}
    HH_uipar = {}
    HH_vipol0 = {}
    CH_Qi = {}
    TG_Qi = {}
    S_jboot = {}
    S_ki = {}
    S_uipar = {}
    S_vipol0 = {}
    HR_PHI = {}
    HS_GAMMA = {}
    HS_Q = {}
    SIM_PHI = {}
    SIM_jboot = {}
    SIM_0vtor0 = {}
    SIM_0upar = {}
    SIM_GAMMA = {}
    SIM_Q = {}
    SIM_PI = {}
    SIM_1upar = {}
    SIM_k = {}
    SIM_K = {}
    SIM_vpol0 = {}
    SIM_1vtor0 = {}
    SIM_upar_by_theta = {}

    Example Usage:
        >>> import matplotlib.pyplot as plt
        >>> from pyrats.data import NEOData
        >>> sim1 = NEOData('example_directory')
        >>> sim1.plot('jboot')
        >>> plt.show()
"""
    # Methods
    def __init__(self, sim_directory):
        """Constructor reads in data from sim_directory and creates new object.

        Every subdirectory of sim_directory will be searched for NEO output
        files, and they will be stored in the self.data dictionary."""

        import numpy
        import os
        import sys

        self.init_data()
        self.master = sim_directory
        self.dirtree = os.walk(sim_directory)
        flag = 0
        for root, dirs, files in self.dirtree:
            if ('out.neo.grid' in files) and ('out.neo.equil' in files) and ('out.neo.theory' in files) and ('out.neo.transport' in files) and ('out.neo.transport_gv' in files):
                flag = 1
                self.toplot.append(root)
                self.set_directory(root)
                self.read_data()
        if flag == 0:
            print "ERROR: NEO files not in directory: " + self.master
            sys.exit()
        self.store_data()

    def init_data(self):
        """Initialize object data."""

        self.master = ""
        self.directory_name = ""
        self.toplot = []
        self.fignum = 1
        self.plotcounter = 1
        self.transport = {}
        self.HH_theory = {}
        self.CH_theory = {}
        self.TG_theory = {}
        self.S_theory = {}
        self.HR_theory = {}
        self.HS_theory = {}
        self.control = {}
        self.n_species = {}
        self.n_energy = {}
        self.n_xi = {}
        self.n_theta = {}
        self.theta_gpoints = {}
        self.n_radial = {}
        self.radial_gpoints = {}
        self.r = {}
        self.dPHI0dr = {}
        self.q = {}
        self.rho_star = {}
        self.rmaj_over_a = {}
        self.omega_rot = {}
        self.omega_rot_deriv = {}
        self.nnorm = {}
        self.Tnorm = {}
        self.n0_over_nnorm = {}
        self.T0_over_Tnorm = {}
        self.a_over_Ln = {}
        self.a_over_LT = {}
        self.inv_tau_self = {}
        self.f = {}
        self.SIM_PHI_by_theta = {}
        self.HH_GAMMA = {}
        self.HH_Qi = {}
        self.HH_Qe = {}
        self.HH_jboot = {}
        self.HH_ki = {}
        self.HH_uipar = {}
        self.HH_vipol0 = {}
        self.CH_Qi = {}
        self.TG_Qi = {}
        self.S_jboot = {}
        self.S_ki = {}
        self.S_uipar = {}
        self.S_vipol0 = {}
        self.HR_PHI = {}
        self.HS_GAMMA = {}
        self.HS_Q = {}
        self.SIM_PHI = {}
        self.SIM_jboot = {}
        self.SIM_0vtor0 = {}
        self.SIM_0upar = {}
        self.SIM_GAMMA = {}
        self.SIM_Q = {}
        self.SIM_PI = {}
        self.SIM_1upar = {}
        self.SIM_k = {}
        self.SIM_K = {}
        self.SIM_vpol0 = {}
        self.SIM_1vtor0 = {}
        self.SIM_upar_by_theta = {}
        self.SIM_GAMMA_gv = {}
        self.SIM_Q_gv = {}
        self.SIM_PI_gv = {}

    def set_directory(self, path):
        """Set the simulation directory."""

        from os.path import expanduser, expandvars
        self.directory_name = expanduser(expandvars(path))
    
    def get_input(self, input_name):
        """Return the specified variable from input.neo.gen.

        input_name  -  requested input

        Ex: get_input("NEO_MODE")
        """

        input_file = file(self.directory_name + 'input.neo.gen', 'r')
        for line in input_file:
            try:
                if line.split()[1] == input_name:
                    return float(line.split()[0])
            except IndexError:
                print "Cannot find specified input parameter: ", input_name
                return 0

    def read_data(self):
        """Read in object data."""
        self.read_grid()
        self.read_equil()
        self.read_theory()
        self.read_transport()
        self.read_transport_gv()

    def store_data(self):
        """Stores data into data dictionary by variable name and directory.

        data can be accessed with two dictionary keys, like so:
        self.data[parameter][directory]."""

        self.transport['PHI'] = NEOOutput(self.SIM_PHI, '(e*PHI1/Tnorm)^2',
                                         'first-order electrostatic potential')
        self.transport['jboot'] = NEOOutput(self.SIM_jboot,
               'jpar*B/(e*nnorm*vnorm*Bunit)', 'first-order bootstrap current')
        self.transport['0vtor0'] = NEOOutput(self.SIM_0vtor0,
                                          '(v_sup_0_sub_phi at theta=0)/vnorm',
      'zeroth-order toroidal flow at the outboard midplane (v_sub_phi = w0*R)')
        self.transport['0upar'] = NEOOutput(self.SIM_0upar,
                                            'upar_sup_0*B/(vnorm*Bunit)',
                              'zeroth-order parallel flow (upar_sup_0=w0*I/B)')
        self.transport['GAMMA'] = NEOOutput(self.SIM_GAMMA,
                    'GAMMA/(nnorm*vnorm)', 'second-order radial particle flux')
        self.transport['Q'] = NEOOutput(self.SIM_Q, 'Q/(nnorm*vnorm*Tnorm)',
                                        'second-order radial energy flux')
        self.transport['PI'] = NEOOutput(self.SIM_PI, 'PI/(nnorm*a*Tnorm)',
                                         'second-order radial momentum flux')
        self.transport['1upar'] = NEOOutput(self.SIM_1upar,
                           'upar*B/(vnorm*Bunit)', 'first-order parallel flow')
        self.transport['k'] = NEOOutput(self.SIM_k, 'k',
                                  'first-order dimensionless flow coefficient')
        self.transport['K'] = NEOOutput(self.SIM_K, 'K/(nnorm*vnorm/Bunit)',
                                    'first-order dimensional flow coefficient')
        self.transport['vpol0'] = NEOOutput(self.SIM_vpol0,
                                            '(v_sub_theta at theta=0)/vnorm',
                          'first-order poloidal flow at the outboard midplane')
        self.transport['1vtor0'] = NEOOutput(self.SIM_1vtor0,
                                             '(v_sub_phi at theta=0)/vnorm',
                          'first-order toroidal flow at the outboard midplane')
        self.control['n_species'] = NEOOutput(self.n_species, 'N_SPECIES',
                                              'number of kinetic species')
        self.control['n_energy'] = NEOOutput(self.n_energy, 'N_ENERGY',
                                             'number of energy polynomials')
        self.control['n_xi'] = NEOOutput(self.n_xi, 'N_XI',
                            'number of xi (cosine of pitch angle) polynomials')
        self.control['n_theta'] = NEOOutput(self.n_theta, 'N_THETA',
                                            'number of theta gridpoints')
        self.control['theta_gpoints'] = NEOOutput(self.theta_gpoints,
                              'theta_sub_j', 'theta gridpoints (j=1..N_THETA)')
        self.control['n_radial'] = NEOOutput(self.n_radial, 'N_RADIAL',
                                             'number of radial gridpoints')
        self.control['radial_gpoints'] = NEOOutput(self.radial_gpoints,
                                                   'r_sub_j/a',
        'radial gridpoints (normalized midplane minor radius) (j=1..N_RADIAL)')
        self.control['r'] = NEOOutput(self.r, 'r/a',
                                      'normalized midplane minor radius')
        self.control['dPHI0dr'] = NEOOutput(self.dPHI0dr,
                                            '(dPHI0dr)(a*e/Tnorm)',
                          'normalized equilibrium-scale radial electric field')
        self.control['q'] = NEOOutput(self.q, 'q', 'safety factor')
        self.control['rho_star'] = NEOOutput(self.rho_star,
                          'rho_star=(c*sqrt(mnorm*Tnorm))*(e*Bunit*a)', 'ratio of the Larmor radius of the normalizing species to the normalizing length scale')
        self.control['rmaj_over_a'] = NEOOutput(self.rmaj_over_a, 'R0/a',
                                 'normailzed flux-surface-center major radius')
        self.control['omega_rot'] = NEOOutput(self.omega_rot, 'w0/(vnorm/a)',
                                       'normailzed toroidal angular frequency')
        self.control['omega_rot_deriv'] = NEOOutput(self.omega_rot_deriv,
                                                    '(dw0/dr)(a^2/vnorm)',
                                          'normalized toroidal rotation shear')
        self.control['nnorm'] = NEOOutput(self.nnorm, 'nnorm',
                             'normalizing equilibrium-scale density (e19/m^3)')
        self.control['Tnorm'] = NEOOutput(self.Tnorm, 'Tnorm',
                             'normalizing equilibrium-scale temperature (keV)')
        self.control['n0_over_nnorm'] = NEOOutput(self.n0_over_nnorm,
                            'n0/nnorm', 'normalized equilibrium-scale density')
        self.control['T0_over_Tnorm'] = NEOOutput(self.T0_over_Tnorm,
                        'T0/Tnorm', 'normalized equilibrium-scale temperature')
        self.control['a_over_Ln'] = NEOOutput(self.a_over_Ln,
                                                   'a/Ln=a(-dln(n0)/dr)',
                                        'normalized equilibrium-scale density')
        self.control['a_over_LT'] = NEOOutput(self.a_over_LT,
                                                   'a/LT=a(-dln(T0)/dr)',
                                    'normalized equilibrium-scale temperature')
        self.control['inv_tau_self'] = NEOOutput(self.inv_tau_self,
                                                 'tau_self^-1/(vnorm/a)',
                                         'normalized self-collision frequency')
        self.HH_theory['GAMMA'] = NEOOutput(self.HH_GAMMA,
                                            'GAMMA/(nnorm*vnorm)',
              'Hinton-Hazeltine second-order radial particle flux (ambipolar)')
        self.HH_theory['Qi'] = NEOOutput(self.HH_Qi, 'Qi/(nnorm*vnorm*Tnorm)',
                      'Hinton-Hazeltine second-order radial energy flux (ion)')
        self.HH_theory['Qe'] = NEOOutput(self.HH_Qe, 'Qe/(nnorm*vnorm*Tnorm)',
                 'Hinton-Hazeltine second-order radial energy flux (electron)')
        self.HH_theory['jboot'] = NEOOutput(self.HH_jboot,
                                            'jpar*B/(e*nnorm*vnorm*Bunit)',
                              'Hinton-Hazeltine first-order bootstrap current')
        self.HH_theory['ki'] = NEOOutput(self.HH_ki, 'ki',
           'Hinton-Hazeltine first-order dimensionless flow coefficient (ion)')
        self.HH_theory['uipar'] = NEOOutput(self.HH_uipar,
                                            'uipar*B/(vnorm*Bunit)',
                            'Hinton-Hazeltine first-order parallel flow (ion)')
        self.HH_theory['vipol0'] = NEOOutput(self.HH_vipol0,
                                             '(vitheta at theta=0)/vnorm',
   'Hinton-Hazeltine first-order poloidal flow at the outboard midplane (ion)')
        self.CH_theory['Qi'] = NEOOutput(self.CH_Qi, 'Qi/(nnorm*vnorm*Tnorm)',
                          'Chang-Hinton second-order radial energy flux (ion)')
        self.TG_theory['Qi'] = NEOOutput(self.TG_Qi, 'Qi/(nnorm*vnorm*Tnorm)',
                               'Taguchi second-order radial energy flux (ion)')
        self.S_theory['jboot'] = NEOOutput(self.S_jboot,
                                           'jpar*B/(e*nnorm*vnorm*Bunit)',
                                        'Sauter first-order bootstrap current')
        self.S_theory['ki'] = NEOOutput(self.S_ki, 'ki',
                     'Sauter first-order dimensionless flow coefficient (ion)')
        self.S_theory['uipar'] = NEOOutput(self.S_uipar,
                                           'uipar*B/(vnorm*Bunit)',
                                      'Sauter first-order parallel flow (ion)')
        self.S_theory['vipol0'] = NEOOutput(self.S_vipol0,
                                            '(vitheta at theta=0)/vnorm',
             'Sauter first-order poloidal flow at the outboard midplane (ion)')
        self.HR_theory['PHI'] = NEOOutput(self.HR_PHI, '(e*PHI1/Tnorm)^2',
                       'Hinton-Rosenbluth first-order electrostatic potential')
        self.HS_theory['GAMMA'] = NEOOutput(self.HS_GAMMA,
                                            'GAMMA/(nnorm*vnorm)',
                           'Hirshman-Sigmar second-order radial particle flux')
        self.HS_theory['Q'] = NEOOutput(self.HS_Q, 'Q/(nnorm*vnorm*Tnorm)',
                             'Hirshman-Sigmar second-order radial energy flux')
        self.transport['GAMMA_gv'] = NEOOutput(self.SIM_GAMMA_gv,
                                               'GAMMA_gv/(nnorm*vnorm)',
                              'Gyro-viscous second-order radial particle flux')
        self.transport['Q_gv'] = NEOOutput(self.SIM_Q_gv,
                                           'Q_gv/(nnorm*vnorm*Tnorm)',
                                'Gyro-viscous second-order radial energy flux')
        self.transport['PI_gv'] = NEOOutput(self.SIM_PI_gv,
                                            'PI_gv/(nnorm*a*Tnorm)',
                              'Gyro-viscous second-order radial momentum flux')

    def read_file(self, fname):
        """Reads NEO data file."""

        import numpy as np

        filename = self.directory_name + '/out.neo.' + fname
        f = np.atleast_2d(np.loadtxt(file(filename)))
        return f

    def read_equil(self):
        """Read out.neo.equil."""

        equil = self.read_file('equil')
        self.r[self.directory_name] = equil[0:self.n_radial[self.directory_name], 0]
        self.dPHI0dr[self.directory_name] = equil[0:self.n_radial[self.directory_name], 1]
        self.q[self.directory_name] = equil[0:self.n_radial[self.directory_name], 2]
        self.rho_star[self.directory_name] = equil[0:self.n_radial[self.directory_name], 3]
        self.rmaj_over_a[self.directory_name] = equil[0:self.n_radial[self.directory_name], 4]
        self.omega_rot[self.directory_name] = equil[0:self.n_radial[self.directory_name], 5]
        self.omega_rot_deriv[self.directory_name] = equil[0:self.n_radial[self.directory_name], 6]
        self.nnorm[self.directory_name] = equil[0:self.n_radial[self.directory_name], 7]
        self.Tnorm[self.directory_name] = equil[0:self.n_radial[self.directory_name], 8]
        self.n0_over_nnorm[self.directory_name] = range(int(self.n_species[self.directory_name]))
        self.T0_over_Tnorm[self.directory_name] = range(int(self.n_species[self.directory_name]))
        self.a_over_Ln[self.directory_name] = range(int(self.n_species[self.directory_name]))
        self.a_over_LT[self.directory_name] = range(int(self.n_species[self.directory_name]))
        self.inv_tau_self[self.directory_name] = range(int(self.n_species[self.directory_name]))
        for a in range(int(self.n_species[self.directory_name])):
            self.n0_over_nnorm[self.directory_name][a] = equil[0:self.n_radial[self.directory_name], 9+5*a]
            self.T0_over_Tnorm[self.directory_name][a] = equil[0:self.n_radial[self.directory_name], 10+5*a]
            self.a_over_Ln[self.directory_name][a] = equil[0:self.n_radial[self.directory_name], 11+5*a]
            self.a_over_LT[self.directory_name][a] = equil[0:self.n_radial[self.directory_name], 12+5*a]
            self.inv_tau_self[self.directory_name][a] = equil[0:self.n_radial[self.directory_name], 13+5*a]

    def read_grid(self):
        """Reads out.neo.grid"""

        import numpy as np

        filename = self.directory_name + '/out.neo.grid'
        grid = np.atleast_1d(np.loadtxt(file(filename)))
        self.n_species[self.directory_name] = grid[0]
        self.n_energy[self.directory_name] = grid[1]
        self.n_xi[self.directory_name] = grid[2]
        self.n_theta[self.directory_name] = grid[3]
        self.theta_gpoints[self.directory_name] = grid[4:4+self.n_theta[self.directory_name]]
        self.n_radial[self.directory_name] = grid[4+self.n_theta[self.directory_name]]
        self.radial_gpoints[self.directory_name] = grid[5+self.n_theta[self.directory_name]:6+self.n_theta[self.directory_name]+self.n_radial[self.directory_name]]

    def read_theory(self):
        """Reads out.neo.theory."""

        theory = self.read_file('theory')
        self.HH_GAMMA[self.directory_name] = theory[0:self.n_radial[self.directory_name], 1]
        self.HH_Qi[self.directory_name] = theory[0:self.n_radial[self.directory_name], 2]
        self.HH_Qe[self.directory_name] = theory[0:self.n_radial[self.directory_name], 3]
        self.HH_jboot[self.directory_name] = theory[0:self.n_radial[self.directory_name], 4]
        self.HH_ki[self.directory_name] = theory[0:self.n_radial[self.directory_name], 5]
        self.HH_uipar[self.directory_name] = theory[0:self.n_radial[self.directory_name], 6]
        self.HH_vipol0[self.directory_name] = theory[0:self.n_radial[self.directory_name], 7]
        self.CH_Qi[self.directory_name] = theory[0:self.n_radial[self.directory_name], 8]
        self.TG_Qi[self.directory_name] = theory[0:self.n_radial[self.directory_name], 9]
        self.S_jboot[self.directory_name] = theory[0:self.n_radial[self.directory_name], 10]
        self.S_ki[self.directory_name] = theory[0:self.n_radial[self.directory_name], 11]
        self.S_uipar[self.directory_name] = theory[0:self.n_radial[self.directory_name], 12]
        self.S_vipol0[self.directory_name] = theory[0:self.n_radial[self.directory_name], 13]
        self.HR_PHI[self.directory_name] = theory[0:self.n_radial[self.directory_name], 14]
        self.HS_GAMMA[self.directory_name] = range(int(self.n_species[self.directory_name]))
        self.HS_Q[self.directory_name] = range(int(self.n_species[self.directory_name]))
        for a in range(int(self.n_species[self.directory_name])):
            self.HS_GAMMA[self.directory_name][a] = theory[0:self.n_radial[self.directory_name], 15+2*a]
            self.HS_Q[self.directory_name][a] = theory[0:self.n_radial[self.directory_name], 16+2*a]

    def read_transport(self):
        """Reads out.neo.transport."""

        transport = self.read_file('transport')
        self.SIM_PHI[self.directory_name] = transport[0:self.n_radial[self.directory_name], 1]
        self.SIM_jboot[self.directory_name] = transport[0:self.n_radial[self.directory_name], 2]
        self.SIM_0vtor0[self.directory_name] = transport[0:self.n_radial[self.directory_name], 3]
        self.SIM_0upar[self.directory_name] = transport[0:self.n_radial[self.directory_name], 4]
        self.SIM_GAMMA[self.directory_name] = range(int(self.n_species[self.directory_name]))
        self.SIM_Q[self.directory_name] = range(int(self.n_species[self.directory_name]))
        self.SIM_PI[self.directory_name] = range(int(self.n_species[self.directory_name]))
        self.SIM_1upar[self.directory_name] = range(int(self.n_species[self.directory_name]))
        self.SIM_k[self.directory_name] = range(int(self.n_species[self.directory_name]))
        self.SIM_K[self.directory_name] = range(int(self.n_species[self.directory_name]))
        self.SIM_vpol0[self.directory_name] = range(int(self.n_species[self.directory_name]))
        self.SIM_1vtor0[self.directory_name] = range(int(self.n_species[self.directory_name]))
        for a in range(int(self.n_species[self.directory_name])):
            self.SIM_GAMMA[self.directory_name][a] = transport[0:self.n_radial[self.directory_name], 5+8*a]
            self.SIM_Q[self.directory_name][a] = transport[0:self.n_radial[self.directory_name], 6+8*a]
            self.SIM_PI[self.directory_name][a] = transport[0:self.n_radial[self.directory_name], 7+8*a]
            self.SIM_1upar[self.directory_name][a] = transport[0:self.n_radial[self.directory_name], 8+8*a]
            self.SIM_k[self.directory_name][a] = transport[0:self.n_radial[self.directory_name], 9+8*a]
            self.SIM_K[self.directory_name][a] = transport[0:self.n_radial[self.directory_name], 10+8*a]
            self.SIM_vpol0[self.directory_name][a] = transport[0:self.n_radial[self.directory_name], 11+8*a]
            self.SIM_1vtor0[self.directory_name][a] = transport[0:self.n_radial[self.directory_name], 12+8*a]

    
    def read_transport_gv(self):
        """Reads out.neo.transport_gv."""
        
        transport_gv = self.read_file('transport_gv')
        self.SIM_GAMMA_gv[self.directory_name] = range(int(self.n_species[self.directory_name]))
        self.SIM_Q_gv[self.directory_name] = range(int(self.n_species[self.directory_name]))
        self.SIM_PI_gv[self.directory_name] = range(int(self.n_species[self.directory_name]))
        for a in range(int(self.n_species[self.directory_name])):
            self.SIM_GAMMA_gv[self.directory_name][a] = transport_gv[0:self.n_radial[self.directory_name], 1+3*a]
            self.SIM_Q_gv[self.directory_name][a] = transport_gv[0:self.n_radial[self.directory_name], 2+3*a]
            self.SIM_PI_gv[self.directory_name][a] = transport_gv[0:self.n_radial[self.directory_name], 3+3*a]

    #-----------------------------------------------#
    # Get data back

    def get_transport(self, var):
        """Get requested transport variable.  Remember that var must be a string.

        Returns a dictionary containing numpy arrays of the requested
        variable paired with their containing directories."""
        import numpy as np
        s = 0
        for k, v in self.transport.iteritems():
            if var == k:
                for k2, v2 in v.data.iteritems():
                    for item in v2:
                        s = s + np.sum(item)
        if s == 0:
            return None
        else:
            return self.transport[var]

    def get_HH_theory(self, var):
        """Get requested Hinton-Hazeltine theory variable."""
 
        s = 0       
        for k, v in self.HH_theory.iteritems():
            if var == k:
                for k2, v2 in v.data.iteritems():
                    for item in v2:
                        s = s + float(item)
        if s == 0:
            return None
        else:
            return self.HH_theory[var]

    def get_CH_theory(self, var):
        """Get requested Chang-Hinton theory variable."""

        s = 0
        for k, v in self.CH_theory.iteritems():
            if var == k:
                for k2, v2 in v.data.iteritems():
                    for item in v2:
                        s = s + float(item)
        if s == 0:
            return None
        else:
            return self.CH_theory[var]

    def get_TG_theory(self, var):
        """Get requested Taguchi theory variable."""

        s = 0
        for k, v in self.TG_theory.iteritems():
            if var == k:
                for k2, v2 in v.data.iteritems():
                    for item in v2:
                        s = s + float(item)
        if s == 0:
            return None
        else:
            return self.TG_theory[var]

    def get_S_theory(self, var):
        """Get requested Sauter theory variable."""

        s = 0
        for k, v in self.S_theory.iteritems():
            if var == k:
                for k2, v2 in v.data.iteritems():
                    for item in v2:
                        s = s + float(item)
        if s == 0:
            return None
        else:
            return self.S_theory[var]

    def get_HR_theory(self, var):
        """Get requested Hinton-Rosenbluth theory variable."""

        s = 0
        for k, v in self.HR_theory.iteritems():
            if var == k:
                for k2, v2 in v.data.iteritems():
                    for item in v2:
                        s = s + float(item)
        if s == 0:
            return None
        else:
            return self.HR_theory[var]

    def get_HS_theory(self, var):
        """Get requested Hirshman-Sigmar theory variable."""

        import numpy as np
        s = 0
        for k, v in self.HS_theory.iteritems():
            if var == k:
                for k2, v2 in v.data.iteritems():
                    for item in v2:
                        s = s + np.sum(item)
        if s == 0:
            return None
        else:
            return self.HS_theory[var]

    def get_control(self, var):
        """Get requested control variable."""

        s = 0
        for k, v in self.control.iteritems():
            if var == k:
                for k2, v2 in v.data.iteritems():
                    for item in v2:
                        s = s + float(item)
        if s == 0:
            return None
        else:
            return self.control[var]

    def print_vars(self):
        """Prints all available variables."""

        for key in sorted(self.data.keys()):
            print key

    #-----------------------------------------------#
    # Plotting routines

    def plot(self, var, n1=2, n2=2, legend=True, verbose=False, cols='bgkcmyrw', styles=['-','--','-.',':']):
        """Plots var as a scatter plot with data from different directories
        coming in different colors."""

        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import os
        import numpy as np
        import sys
        #Adjust subplot display parameters
        mpl.rcParams['font.size'] = 10.0
        mpl.rcParams['figure.subplot.wspace'] = .99
        mpl.rcParams['figure.subplot.hspace'] = .3
        mpl.rcParams['figure.subplot.right'] = .8

        #Check to see if requested variable is actually available
        flag = 0
        for k in self.transport.iterkeys():
            if var == k:
                flag = 1
        if not flag:
            print "ERROR: ", var, " is not a valid parameter.  Type -options for help."
            sys.exit()

        #Add related variables to todo list
        varlist = []
        varlist.append(var)
        if var == 'Q':
            varlist.append('Qi')
            varlist.append('Qe')
        if var == 'k':
            varlist.append('ki')
        if var == '1upar':
            varlist.append('uipar')
        if var == 'vpol0':
            varlist.append('vipol0')

        #Plotcounter keeps track of where on the current figure the current plot
        #is
        if self.plotcounter > (n1 * n2):
            self.plotcounter = 1
            self.fignum = self.fignum + 1
        fig = plt.figure(self.fignum)
        ax = fig.add_subplot(n1, n2, self.plotcounter)
        ax.set_xlabel('r/a')
        ax.set_ylabel(self.get_transport(var).units)
        #verbose sets title to be more descriptive
        if verbose:
            ax.set_title(self.get_transport(var).descriptor+' vs. r/a')
        else:
            ax.set_title(var + ' vs. r/a')
        for var in varlist:
            tempr = []
            #Obtain all radius values
            for key in self.toplot:
                tempr.append(self.get_control('r').data[key])
            #Store all values in a 1-D numpy array
            tempr = np.array(self.split(tempr)).flatten()
            #Sort them by increasing radius
            ind = tempr.argsort()
            tempr = tempr[ind]
            print tempr

            #If the requested variable is available
            if self.get_transport(var) != None:
                transport = []
                #find it
                for key in self.toplot:
                    transport.append(self.get_transport(var).data[key])
                #Check to see if it has one or two dimensions
                try:
                    transport[0][0][0]
                    #If we reach this point, that means it has two dimensions,
                    #so we have to sort that out properly
                    for y in range(len(transport[0])):
                        plot = []
                        for x in range(len(transport)):
                            plot.append(transport[x][y])
                        plot = np.array(self.split(plot)).flatten()[ind]
                        #Create plot of sorted data
                        ax.plot(tempr,plot,c=cols[0],ls=styles[y],label='Sim spe '+str(y))
                except IndexError:
                    #If we reach this point, that means that there is only one
                    #dimension.  All we have to do is create the array, sort it,
                    #and plot it.
                    transport = np.array(self.split(transport)).flatten()[ind]
                    print transport
                    ax.plot(tempr, transport, c=cols[0], label='Sim')

            if self.get_HH_theory(var) != None:
                #HH_theory always has only one dimension
                HH = []
                for key in self.toplot:
                    HH.append(self.get_HH_theory(var).data[key])
                print var
                print np.array(self.split(HH)).flatten()
                print ind
                HH = np.array(self.split(HH)).flatten()[ind]
                ax.plot(tempr, HH, c=cols[1], label='HH ' + var)

            if self.get_CH_theory(var) != None:
                #So does CH_theory
                CH = []
                for key in self.toplot:
                    CH.append(self.get_CH_theory(var).data[key])
                CH = np.array(self.split(CH)).flatten()[ind]
                ax.plot(tempr, CH, c=cols[2], label='CH ' + var)

            if self.get_TG_theory(var) != None:
                #And TG_theory
                TG = []
                for key in self.toplot:
                    TG.append(self.get_TG_theory(var).data[key])
                TG = np.array(self.split(TG)).flatten()[ind]
                ax.plot(tempr, TG, c=cols[3], label='TG ' + var)

            if self.get_S_theory(var) != None:
                #etc
                S = []
                for key in self.toplot:
                    S.append(self.get_S_theory(var).data[key])
                S = np.array(self.split(S)).flatten()[ind]
                ax.plot(tempr, S, c=cols[4], label='S ' + var)

            if self.get_HR_theory(var) != None:
                HR = []
                for key in self.toplot:
                    HR.append(self.get_HR_theory(var).data[key])
                HR = np.array(self.split(HR)).flatten()[ind]
                ax.plot(tempr, HR, c=cols[5], label='HR ' + var)

            if self.get_HS_theory(var) != None:
                #HS_theory always has two dimensions, so we have to sort that
                #out the same way we sorted out transport
                HS = []
                for key in self.toplot:
                    HS.append(self.get_HS_theory(var).data[key])
                for y in range(len(HS[0])):
                    plot = []
                    for x in range(len(HS)):
                        plot.append(HS[x][y])
                    plot = np.array(self.split(plot)).flatten()[ind]
                    ax.plot(tempr,plot,c=cols[6],ls=styles[y],label='HS '+var+' spe '+str(y))
        #If the legend has been requested, then we create it
        if legend == True:
            ax.legend(loc=2,bbox_to_anchor=(1,1))
        #Move on to the next plot
        self.plotcounter = self.plotcounter + 1

    #----------------------------------------------------------#
    # Misc methods

    def split(self, array):
        """split splits an array which may be made up of elements with
        multiple entries into an array with one entry per element."""

        import numpy as np

        t = []
        for s in array:
            if len(s) > 1:
                for a in s:
                    t.append(np.array([a]))
            else:
                t.append(s)
        return t

#==============================================================================
#==============================================================================

class GYROData:
    """A class of GYRO output data.

    Data:

    directory_name = ""
    profile = {}
    geometry = {}
    t = {}
    freq = {}
    diff = {}
    diff_i = {}
    diff_n = {}
    gbflux = []
    gbflux_i = []
    gbflux_n = {}
    moment_u = []
    moment_n = []
    moment_e = []
    moment_v = []
    moment_zero = []
    flux_velocity = {}
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
        path = '~/gacode/shared/python/pyrats'
        sys.path.append(expanduser(expandvars(path)))

    def init_data(self):
        """Initialize object data."""

        self.directory_name = ""
        self.profile = {}
        self.geometry = {}
        self.t = {}
        self.freq = {}
        self.diff = {}
        self.diff_i = {}
        self.diff_n = {}
        self.gbflux = []
        self.gbflux_i = []
        self.gbflux_n = {}
        self.moment_u = []
        self.moment_n = []
        self.moment_e = []
        self.moment_v = []
        self.moment_zero = []
        self.flux_velocity = {}
        self.k_perp_squared = []

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
        #self.read_diff()
        #self.read_diff_i()
        #self.read_diff_n()
        self.read_gbflux()
        #self.read_gbflux_i()
        #self.read_gbflux_n()
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
        import fileupload
        import sys
        import os

        fname = 'out.gyro.' + fname
        filename = self.directory_name + '/' + fname
        if fname in os.listdir(self.directory_name):
            f = fileupload.loadtxt(filename, dSize)
            return f
        else:
            print "ERROR: File " + fname + " does not exist."
            return 0


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
        if geometry != 0:
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

    def read_diff(self):
        """Reads in diff data.  Output is dictionary of numpy arrays with
        dimensions: n_kinetic x n_field x n_time"""

        import numpy as np
    
        diff = self.read_file('diff', 12)
        if diff != 0:
            temp = diff.reshape( (self.t['n_time'], self.profile['n_kinetic'], self.profile['n_field'], 2), order='F')
            self.diff['D_sigma/chi_{GB}'] = temp[:, :, :, 0]
            self.diff['chi_sigma/chi_{GB}'] = temp[:, :, :, 1]

    def read_diff_i(self):
        """Reads in diff_i.  Output is dictionary of numpy arrays with
        dimensions: n_kinetic x n_field x n_x x n_time"""

        import numpy as np

        diff_i = self.read_file('diff_i', 12)
        if diff_i != 0:
            temp = diff_i.reshape( (self.t['n_time'], self.profile['n_x'], 2, self.profile['n_field'], self.profile['n_kinetic']), order='F')
            self.diff_i['D_sigma(r)/chi_{GB}'] = temp[:, :, 0, :, :]
            self.diff_i['chi_sigma(r)/chi_{GB}'] = temp[:, :, 1, :, :]

    def read_diff_n(self):
        """Reads in diff_n.  Output is dictionary of numpy arrays with
        dimensions: n_kinetic x n_field x n_n x n_time"""

        import numpy as np

        diff_n = self.read_file('diff_n', 12)
        if diff_n != 0:
            temp = diff_n.reshape( (self.t['n_time'], self.profile['n_n'], 2, self.profile['n_field'], self.profile['n_kinetic']), order='F')
            self.diff_n['D_{sigma,n}/chi_{GB}'] = temp[:, :, 0, :, :]
            self.diff_n['chi_{sigma,n}/chi_{GB}'] = temp[:, :, 1, :, :]

    def read_gbflux(self):
        """Reads in gbflux data.  Output is numpy array with dimensions:
        n_kinetic x n_field x 4 x n_time"""

        import numpy as np

        gbflux = self.read_file('gbflux', 12)
        if gbflux != 0:
            self.gbflux = gbflux.reshape( (self.t['n_time'], self.profile['n_kinetic'], self.profile['n_field'], 4), order='F')      

    def read_gbflux_i(self):
        """Reads in gbflux_i data.  Output is numpy array with dimensions:
        n_kinetic x n_field x 4 x n_x x n_time"""

        import numpy as np

        gbflux_i = self.read_file('gbflux_i.conv', 12)
        if gbflux_i != 0:
            self.gbflux_i = gbflux_i.reshape( (self.profile['n_kinetic'], self.profile['n_field'], 4, self.profile['n_x'], self.t['n_time']), order='F')

    def read_gbflux_n(self):
        """Reads gbflux_n data.  Output is numpy array with dimensions:
        n_kinetic x n_field x 4 x n_n x n_time"""

        import numpy as np

        gbflux_n = self.read_file('gbflux_n', 12)
        if gbflux_n != 0:
            temp = gbflux_n.reshape( (self.profile['n_kinetic'], self.profile['n_field'], 4, self.profile['n_n'], self.t['n_time']), order='F')
            self.gbflux_n['GAMMA_sigma(r)/GAMMA_{GB}'] = temp[:, :, 0, :, :]
            self.gbflux_n['Q_{0,sigma}(r)/Q_{GB}'] = temp[:, :, 1, :, :]
            self.gbflux_n['PI_sigma(r)/PI_{GB}'] = temp[:, :, 2, :, :]
            self.gbflux_n['S^tur_{W,sigma}(r)/S_{GB}'] = temp[:, :, 3, :, :]

    def read_moment_u(self):
        """Reads in moment_u data.  Output is numpy array with dimensions:
        2 x n_theta_plot x n_x x n_field x n_n x n_time"""

        import numpy as np
        import time
        import fileupload

        starttime = time.time()
        moment_u = self.read_file('moment_u', 12)
        if moment_u != 0:
            self.moment_u = moment_u.reshape( (2, self.profile['n_theta_plot'], self.profile['n_x'], self.profile['n_field'], self.profile['n_n'], self.t['n_time']), order='F')
            endtime = time.time()
            print "Time: " + str(endtime - starttime)

    def read_moment_n(self):
        """Reads in moment_n data.  Output is numpy array with dimensions:
        2 x n_theta_plot x n_x x n_kinetic x n_n x n_time"""

        import numpy as np
        import time
        import fileupload

        starttime = time.time()
        moment_n = self.read_file('moment_n', 12)
        if moment_n != 0:
            self.moment_n = moment_n.reshape( (2, self.profile['n_theta_plot'], self.profile['n_x'], self.profile['n_kinetic'], self.profile['n_n'], self.t['n_time']), order='F')
            endtime = time.time()
            print "Time: " + str(endtime - starttime)

    def read_moment_e(self):
        """Reads in moment_e data.  Output is numpy array with dimensions:
        2 x n_theta_plot x n_x x n_kinetic x n_n x n_time"""

        import numpy as np
        import time
        import fileupload

        starttime = time.time()
        moment_e = self.read_file('moment_e', 12)
        if moment_e != 0:
            self.moment_e = moment_e.reshape( (2, self.profile['n_theta_plot'], self.profile['n_x'], self.profile['n_kinetic'], self.profile['n_n'], self.t['n_time']), order='F')
            endtime = time.time()
            print "Time: " + str(endtime - starttime)

    def read_moment_v(self):
        """Reads in moment_v data.  Output is numpy array with dimensions:
        2 x n_theta_plot x n_x x n_kinetic x n_n x n_time"""

        import numpy as np
        import time
        import fileupload

        starttime = time.time()
        moment_v = self.read_file('moment_v', 12)
        if moment_v != 0:
            midtime = time.time()
            self.moment_v = moment_v.reshape( (2, self.profile['n_theta_plot'], self.profile['n_x'], self.profile['n_kinetic'], self.profile['n_n'], self.t['n_time']), order='F')
            endtime = time.time()
            print "Time to load data: " + str(midtime - starttime)
            print "Time to reshape data: " + str(endtime - midtime)
            print "Total time: " + str(endtime - starttime)

    def read_moment_zero(self):
        """Reads in moment_zero data.  Output is numpy array with dimensions:
        n_x x n_kinetic x n_moment x n_time"""

        import numpy as np
        import time

        starttime = time.time()
        moment_zero = self.read_file('moment_zero', 12)
        if moment_zero != 0:
            self.moment_zero = moment_zero.reshape( (self.profile['n_x'], self.profile['n_kinetic'], self.profile['n_moment'], self.t['n_time']), order='F')
            endtime = time.time()
            print "Time: " + str(endtime - starttime)

    def read_flux_velocity(self):
        """Reads in flux_velocity data.  Output is numpy array with dimensions:
        n_energy x n_lambda x n_kinetic x n_field x 2 x n_n x n_time"""

        import numpy as np

        flux_velocity = self.read_file('flux_velocity', 12)
        if flux_velocity != 0:
            temp = flux_velocity.reshape( (self.t['n_time'], self.profile['n_n'], 2, self.profile['n_field'], self.profile['n_kinetic'], self.profile['n_lambda'], self.profile['n_energy']), order='F')
            self.flux_velocity['GAMMA_{sigma, n}(epsilon, lambda)'] = temp[:, :, 0, :, :, :, :]
            self.flux_velocity['Q_{sigma, n}(epsilon, lambda)'] = temp[:, :, 1, :, :, :, :]

    def read_k_perp_squared(self):
        """Reads in k_perp_squared data.  Output is numpy array with dimensions:
        n_n x n_time"""

        import numpy as np

        k_perp_squared = self.read_file('k_perp_squared', 12)
        if k_perp_squared != 0:
            self.k_perp_squared = k_perp_squared.reshape( (self.t['n_time'], self.profile['n_n']), order='F')

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
        for i in range(int(self.profile['n_spec'])):
            temp1.append(self.gbflux[i, :, 0, :]*np.mean(self.profile['dlnndr_s'], axis=1)[i]/np.mean(self.profile['den_s'], axis=1)[i])
            temp2.append(self.gbflux[i, :, 1, :]*np.mean(self.profile['dlntdr_s'], axis=1)[i]/np.mean(self.profile['tem_s'], axis=1)[i])
        self.diff['D_sigma/chi_{GB}'] = np.array(temp1)
        self.diff['chi_sigma/chi_{GB}'] = np.array(temp2)

    def make_diff_i(self):
        """Makes diff_i.  Output is dictionary of numpy arrays with
        dimensions: n_kinetic x n_field x n_x x n_time"""

        import numpy as np

        if self.gbflux_i == []:
            self.read_gbflux_i()
        temp1 = []
        temp2 = []
        for i in range(int(self.profile['n_spec'])):
            temp1.append(self.gbflux_i[i, :, 0, :, :]*np.mean(self.profile['dlnndr_s'], axis=1)[i]/np.mean(self.profile['den_s'], axis=1)[i])
            temp2.append(self.gbflux[i, :, 1, :, :]*np.mean(self.profile['dlntdr_s'], axis=1)[i]/np.mean(self.profile['tem_s'], axis=1)[i])
        self.diff_i['D_sigma(r)/chi_{GB}'] = np.array(temp1)
        self.diff_i['chi_sigma(r)/chi_{GB}'] = np.array(temp2)

    def plot(self, x, y, dim=(1,1)):
        """Creates matplotlib plot of requested data."""

        import matplotlib.pyplot as plt
        import numpy as np

        fig = plt.figure(self.fignum)
        self.fignum = self.fignum + 1
        ax = fig.add_subplot(dim[0], dim[1], self.plotcounter)
        ax.plot(x, y)
