class profiles_genData:
    """A class of profiles_gen data.

     Data:

"""


    #Methods
    def __init__(self, directory = '.'):
        """Constructor reads in data from directory and creates new object."""
        
        self.set_directory(directory)
        self.init_data()
        self.store_data()

    def init_data(self):
        """Initialize object data."""

        self.data = []
        self.n_exp = 0
        self.hlen = 0
        self.ar = []
        self.br = []
        self.az = []
        self.bz = []

    def set_directory(self, directory):
        """Set the directory which contains input.profiles."""

        from os.path import expanduser, expandvars
        path = directory
        self.directory_name = expanduser(expandvars(path))

    def read_data(self):
        """Read in object data from input.profiles."""

        import numpy as np

        elements = {}
        temp = []
        raw_data = open(self.directory_name + '/input.profiles', 'r').readlines()
        while raw_data[self.hlen].strip()[0].isdigit() == False:
            self.hlen = self.hlen + 1
            if raw_data[self.hlen].strip()[0:6] == 'N_EXP=':
                self.n_exp = int(raw_data[self.hlen].strip()[6:])
        for line in range(self.hlen, len(raw_data)):
            if raw_data[line].strip()[0].isdigit():
                temp.append(raw_data[line].split())
        self.hlen = self.hlen - 1
        data = np.array(temp)

        keywords = []
        for count in range(5):
            keywords = raw_data[(self.n_exp + 2)*count + self.hlen].replace('(kW/eV)', '(kW/eV) ').split()
            column = 0
            for key in keywords:
                elements[key] = data[0:self.n_exp, column]
                column = column + 1
            data = data[self.n_exp:, :]
        for count in range(3):
            keys = raw_data[(self.n_exp + 2) * count + 5*(self.n_exp + 2) + self.hlen].split()
            keywords = []
            for x in range(5):
                keywords.append(keys[2 * x])
            column = 0
            for key in keywords:
                elements[key] = data[0:self.n_exp, column]
                column = column + 1
            data = data[self.n_exp:, :]
        return elements

    def read_fourier(self):
        """Read in data from input.profiles.geo."""

        temp = []
        raw_data = open(self.directory_name + '/input.profiles.geo', 'r').readlines()
        for line in raw_data:
            temp.append(line.split()[0])
        count = int(temp[0]) + 1
        x = 1
        for x in range(self.n_exp):
            arn = []
            brn = []
            azn = []
            bzn = []
            for y in range(count):
                arn.append(temp[4*(count*x + y) + 1])
                brn.append(temp[4*(count*x + y) + 2])
                azn.append(temp[4*(count*x + y) + 3])
                bzn.append(temp[4*(count*x + y) + 4])
            self.ar.append(arn)
            self.br.append(brn)
            self.az.append(azn)
            self.bz.append(bzn)

    def store_data(self):
        """Reads data and renames it appropriately."""

        self.data = self.read_data()
        self.data['rho'] = self.data['#rho(-)']
        self.data['rho (-)'] = self.data['rho']
        self.data.pop('#rho(-)')
        self.data['rmin'] = self.data['rmin(m)']
        self.data['rmin (m)'] = self.data['rmin(m)']
        self.data.pop('rmin(m)')
        self.data['rmaj'] = self.data['rmaj(m)']
        self.data['rmaj (m)'] = self.data['rmaj(m)']
        self.data.pop('rmaj(m)')
        self.data['q'] = self.data['q(-)']
        self.data['q (-)'] = self.data['q(-)']
        self.data.pop('q(-)')
        self.data['kappa'] = self.data['kappa(-)']
        self.data['kappa (-)'] = self.data['kappa(-)']
        self.data.pop('kappa(-)')
        self.data['delta'] = self.data['#delta(-)']
        self.data['delta (-)'] = self.data['delta']
        self.data.pop('#delta(-)')
        self.data['Te'] = self.data['Te(keV)']
        self.data['Te (keV)'] = self.data['Te(keV)']
        self.data.pop('Te(keV)')
        self.data['ne'] = self.data['ne(10^19/m^3)']
        self.data['ne (10^19/m^3)'] = self.data['ne(10^19/m^3)']
        self.data.pop('ne(10^19/m^3)')
        self.data['z_eff'] = self.data['z_eff(-)']
        self.data['z_eff (-)'] = self.data['z_eff(-)']
        self.data.pop('z_eff(-)')
        self.data['omega0'] = self.data['omega0(1/s)']
        self.data['omega0 (1/s)'] = self.data['omega0(1/s)']
        self.data.pop('omega0(1/s)')
        self.data['flow_mom'] = self.data['#flow_mom(N-m)']
        self.data['flow_mom (N-m)'] = self.data['flow_mom']
        self.data.pop('#flow_mom(N-m)')
        self.data['pow_e'] = self.data['pow_e(MW)']
        self.data['pow_e (MW)'] = self.data['pow_e(MW)']
        self.data.pop('pow_e(MW)')
        self.data['pow_i'] = self.data['pow_i(MW)']
        self.data['pow_i (MW)'] = self.data['pow_i(MW)']
        self.data.pop('pow_i(MW)')
        self.data['pow_ei'] = self.data['pow_ei(MW)']
        self.data['pow_ei (MW)'] = self.data['pow_ei(MW)']
        self.data.pop('pow_ei(MW)')
        self.data['zeta'] = self.data['zeta(-)']
        self.data['zeta (-)'] = self.data['zeta(-)']
        self.data.pop('zeta(-)')
        self.data['flow_beam'] = self.data['#flow_beam(kW/eV)']
        self.data['flow_beam (kW/eV)'] = self.data['flow_beam']
        self.data.pop('#flow_beam(kW/eV)')
        self.data['flow_wall'] = self.data['flow_wall(kW/eV)']
        self.data['flow_wall (kW/eV)'] = self.data['flow_wall(kW/eV)']
        self.data.pop('flow_wall(kW/eV)')
        self.data['zmag'] = self.data['zmag(m)']
        self.data['zmag (m)'] = self.data['zmag(m)']
        self.data.pop('zmag(m)')
        self.data['ptot'] = self.data['ptot(Pa)']
        self.data['ptot (Pa)'] = self.data['ptot(Pa)']
        self.data.pop('ptot(Pa)')
        self.data['ni_1'] = self.data['#ni_1(10^19/m^3)']
        self.data['ni_1 (10^19/m^3)'] = self.data['ni_1']
        self.data.pop('#ni_1(10^19/m^3)')
        self.data['ni_2'] = self.data['ni_2(10^19/m^3)']
        self.data['ni_2 (10^19/m^3)'] = self.data['ni_2(10^19/m^3)']
        self.data.pop('ni_2(10^19/m^3)')
        self.data['ni_3'] = self.data['ni_3(10^19/m^3)']
        self.data['ni_3 (10^19/m^3)'] = self.data['ni_3(10^19/m^3)']
        self.data.pop('ni_3(10^19/m^3)')
        self.data['ni_4'] = self.data['ni_4(10^19/m^3)']
        self.data['ni_4 (10^19/m^3)'] = self.data['ni_4(10^19/m^3)']
        self.data.pop('ni_4(10^19/m^3)')
        self.data['ni_5'] = self.data['ni_5(10^19/m^3)']
        self.data['ni_5 (10^19/m^3)'] = self.data['ni_5(10^19/m^3)'] 
        self.data.pop('ni_5(10^19/m^3)')
        self.data['Ti_1'] = self.data['#Ti_1']
        self.data['Ti_1 (keV)'] = self.data['Ti_1']
        self.data.pop('#Ti_1')
        self.data['Ti_2 (keV)'] = self.data['Ti_2']
        self.data['Ti_3 (keV)'] = self.data['Ti_3']
        self.data['Ti_4 (keV)'] = self.data['Ti_4']
        self.data['Ti_5 (keV)'] = self.data['Ti_5']
        self.data['vtor_1'] = self.data['#vtor_1']
        self.data['vtor_1 (m/s)'] = self.data['vtor_1']
        self.data.pop('#vtor_1')
        self.data['vtor_2 (m/s)'] = self.data['vtor_2']        
        self.data['vtor_3 (m/s)'] = self.data['vtor_3']
        self.data['vtor_4 (m/s)'] = self.data['vtor_4']  
        self.data['vtor_5 (m/s)'] = self.data['vtor_5'] 
        self.data['vpol_1'] = self.data['#vpol_1']
        self.data['vpol_1 (m/s)'] = self.data['vpol_1']
        self.data.pop('#vpol_1')
        self.data['vpol_2 (m/s)'] = self.data['vpol_2']
        self.data['vpol_3 (m/s)'] = self.data['vpol_3']
        self.data['vpol_4 (m/s)'] = self.data['vpol_4']
        self.data['vpol_5 (m/s)'] = self.data['vpol_5']

    def compute_mtypeeq(self, r):

        """Uses input data to compute the Miller-type equilibrium."""
        
        import math

        R = []
        Z = []
        tot = []
        theta = 0
        dtheta = 0.01
        x = self.match(r, self.data['rho'])
        while theta < 2 * math.pi:
            a = float(self.data['rmaj'][x]) + float(self.data['rmin'][x]) * math.cos(theta + math.asin(float(self.data['delta'][x])) * math.sin(theta))
            R.append(a)
            b = float(self.data['zmag'][x]) + float(self.data['kappa'][x]) * float(self.data['rmin'][x]) * math.sin(theta + float(self.data['zeta'][x]) * math.sin(2 * theta))
            Z.append(b)
            theta = theta + dtheta
        tot.append(R)
        tot.append(Z)
        tot.append(self.data['rmaj'][x])
        tot.append(self.data['zmag'][x])
        return tot

    def compute_fouriereq(self, r):

        """Uses input data to compute the general Fourier-series equilibrium."""

        import math
        
        self.read_fourier()
        R = []
        Z = []
        tot = []      
        theta = 0
        dtheta = 0.01
        x = self.match(r, self.data['rho'])
        while theta < 2 * math.pi:
            tempr = []
            tempz = []
            for n in range(1, len(self.ar[0])):
                tempr.append(float(self.ar[x][n])*math.cos(n*theta) + float(self.br[x][n])*math.sin(n*theta))
                tempz.append(float(self.az[x][n])*math.cos(n*theta) + float(self.bz[x][n])*math.sin(n*theta))
            a = float(self.ar[x][0])/2 + sum(tempr)
            b = float(self.az[x][0])/2 + sum(tempz)
            R.append(a)
            Z.append(b)
            theta = theta + dtheta
        tot.append(R)
        tot.append(Z)
        tot.append(self.data['rmaj'][x])
        tot.append(self.data['zmag'][x])
        return tot

    #-------------------------------------------- #
    # Get data back
    def get(self, var):
        """Return requested variable."""

        return self.data[var]

    #-------------------------------------------- #
    # Useful Functions
    def match(self, val, vec):
        """Return closest match to input in a list of values."""
        import math

        temp = []
        for n in range(len(vec)):
            t = math.fabs(val - float(vec[n])), n
            temp.append(t)
        return sorted(temp)[0][1]
