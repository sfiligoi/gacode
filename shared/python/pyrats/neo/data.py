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
        files, and they will be stored in either the self.transport or one of
        the self.theory dictionaries."""

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
