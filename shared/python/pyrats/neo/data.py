class NEOData:
    """A class for management of NEO output data.

    Data (things that actually get used):

    master = ""
    dirname = ""
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
    T0_over_Tnorm = {}
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
        >>> from pyrats.neo.data import NEOData
        >>> sim1 = NEOData('example_directory')
    """
    
    #---------------------------------------------------------------------------#
    # Methods

    def __init__(self, sim_directory):
        """Constructor reads in data from sim_directory and creates new object.
        """

        self.init_data()
        self.set_directory(sim_directory)
        self.read_grid()
        self.read_equil()
        self.read_theory()
        self.read_transport()
        self.read_transport_gv()

    #---------------------------------------------------------------------------#

    def init_data(self):
        """Initialize object data."""

        self.master      = ""
        self.dirname     = ""
        self.toplot      = []
        self.fignum      = 1
        self.plotcounter = 1
        self.transport   = {}
        self.HH_theory   = {}
        self.CH_theory   = {}
        self.TG_theory   = {}
        self.S_theory    = {}
        self.HR_theory   = {}
        self.HS_theory   = {}
        self.control     = {}

        # 

        self.n_species = []
        self.n_energy = []
        self.n_xi = []
        self.n_theta = []
        self.theta_gpoints = []
        self.n_radial = []
        self.radial_gpoints = []
        self.r = []
        self.dPHI0dr = []
        self.q = []
        self.rho_star = []
        self.rmaj_over_a = []
        self.omega_rot = []
        self.omega_rot_deriv = []
        self.nnorm = []
        self.Tnorm = []
        self.vnorm_over_a = []
        self.n0_over_nnorm = []
        self.T0_over_Tnorm = []
        self.a_over_Ln = []
        self.a_over_LT = []
        self.inv_tau_self = []
        self.f = []
        self.SIM_PHI_by_theta = []
        self.HH_GAMMA = []
        self.HH_Qi = []
        self.HH_Qe = []
        self.HH_jboot = []
        self.HH_ki = []
        self.HH_uipar = []
        self.HH_vipol0 = []
        self.CH_Qi = []
        self.TG_Qi = []
        self.S_jboot = []
        self.S_ki = []
        self.S_uipar = []
        self.S_vipol0 = []
        self.HR_PHI = []
        self.HS_GAMMA = []
        self.HS_Q = []
        self.SIM_PHI = []
        self.SIM_jboot = []
        self.SIM_0vtor0 = []
        self.SIM_0upar = []
        self.SIM_GAMMA = []
        self.SIM_Q = []
        self.SIM_PI = []
        self.SIM_1upar = []
        self.SIM_k = []
        self.SIM_K = []
        self.SIM_vpol0 = []
        self.SIM_1vtor0 = []
        self.SIM_upar_by_theta = []
        self.SIM_GAMMA_gv = []
        self.SIM_Q_gv = []
        self.SIM_PI_gv = []

    #---------------------------------------------------------------------------#

    def set_directory(self, path):
        """Set the simulation directory."""

        from os.path import expanduser, expandvars
        self.dirname = expanduser(expandvars(path))
    
    #---------------------------------------------------------------------------#

    def read_grid(self):
        """Reads out.neo.grid"""

        import numpy as np
        import sys

        import sys
        import numpy as np
 
        try:
            data = np.loadtxt(self.dirname+'/out.neo.grid')
        except:
            print "ERROR (NEOData): Fatal error!  Missing out.neo.grid."
            sys.exit()
                
        self.n_species.append(grid[0])
        self.n_energy.append(grid[1])
        self.n_xi.append(grid[2])
        self.n_theta.append(grid[3])
        self.theta_gpoints.append(grid[4:4+self.n_theta[-1]])
        self.n_radial.append(grid[4+self.n_theta[-1]])
        self.radial_gpoints.append(grid[5+self.n_theta[-1]:6+self.n_theta[-1]+self.n_radial[-1]])

    def store_data(self):
        """Stores data into data dictionaries by variable name and directory.

        data can be accessed with two dictionary keys, like so:
        self.transport[parameter][directory]."""

        import numpy as np

        self.transport['PHI'] = NEOOutput(np.array(self.SIM_PHI), '(e*PHI1/Tnorm)^2',
                                         'first-order electrostatic potential')
        self.transport['jboot'] = NEOOutput(np.array(self.SIM_jboot),
               'j_par*B/(e*nnorm*vnorm*B_unit)', 'first-order bootstrap current')
        self.transport['0vtor0'] = NEOOutput(np.array(self.SIM_0vtor0),
                                          '(v^{(0)}_phi (theta=0)/vnorm',
      'zeroth-order toroidal flow at the outboard midplane (v^(0)_phi = w0*R)')
        self.transport['0upar'] = NEOOutput(np.array(self.SIM_0upar),
                                            'upar^(0)*B/(vnorm*B_unit)',
                              'zeroth-order parallel flow (upar^(0)_parallel=w0*I/B)')
        self.transport['GAMMA'] = NEOOutput(np.array(self.SIM_GAMMA),
                    'GAMMA/(nnorm*vnorm)', 'second-order radial particle flux')
        self.transport['Q'] = NEOOutput(np.array(self.SIM_Q), 'Q/(nnorm*vnorm*Tnorm)',
                                        'second-order radial energy flux')
        self.transport['PI'] = NEOOutput(np.array(self.SIM_PI), 'PI/(nnorm*a*Tnorm)',
                                         'second-order radial momentum flux')
        self.transport['1upar'] = NEOOutput(np.array(self.SIM_1upar),
                           'upar*B/(vnorm*Bunit)', 'first-order parallel flow')
        self.transport['k'] = NEOOutput(np.array(self.SIM_k), 'k',
                                  'first-order dimensionless flow coefficient')
        self.transport['K'] = NEOOutput(np.array(self.SIM_K), 'K/(nnorm*vnorm/Bunit)',
                                    'first-order dimensional flow coefficient')
        self.transport['vpol0'] = NEOOutput(np.array(self.SIM_vpol0),
                                            '(v_sub_theta at theta=0)/vnorm',
                          'first-order poloidal flow at the outboard midplane')
        self.transport['1vtor0'] = NEOOutput(np.array(self.SIM_1vtor0),
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
        self.control['r'] = NEOOutput(np.array(self.r), 'r/a',
                                      'normalized midplane minor radius')
        self.control['dPHI0dr'] = NEOOutput(np.array(self.dPHI0dr),
                                            '(dPHI0dr)(a*e/Tnorm)',
                          'normalized equilibrium-scale radial electric field')
        self.control['q'] = NEOOutput(np.array(self.q), 'q', 'safety factor')
        self.control['rho_star'] = NEOOutput(np.array(self.rho_star),
                          'rho_star=(c*sqrt(mnorm*Tnorm))*(e*Bunit*a)', 'ratio of the Larmor radius of the normalizing species to the normalizing length scale')
        self.control['rmaj_over_a'] = NEOOutput(np.array(self.rmaj_over_a), 'R0/a',
                                 'normailzed flux-surface-center major radius')
        self.control['omega_rot'] = NEOOutput(np.array(self.omega_rot), 'w0/(vnorm/a)',
                                       'normailzed toroidal angular frequency')
        self.control['omega_rot_deriv'] = NEOOutput(np.array(self.omega_rot_deriv),
                                                    '(dw0/dr)(a^2/vnorm)',
                                          'normalized toroidal rotation shear')
        self.control['nnorm'] = NEOOutput(np.array(self.nnorm), '(e19/m^3)',
                             'normalizing equilibrium-scale density (e19/m^3)')
        self.control['Tnorm'] = NEOOutput(np.array(self.Tnorm), 'keV',
                             'normalizing equilibrium-scale temperature (keV)')
        self.control['vnorm_over_a'] = NEOOutput(np.array(self.vnorm_over_a), '1/s',
                                                 'ratio of the normalizing equilibrium-scale thermal speed to the normalizing length scale (1/s) ')
        self.control['n0_over_nnorm'] = NEOOutput(np.array(self.n0_over_nnorm),
                            'unitless', 'normalized equilibrium-scale density')
        self.control['T0_over_Tnorm'] = NEOOutput(np.array(self.T0_over_Tnorm),
                        'unitless', 'normalized equilibrium-scale temperature')
        self.control['a_over_Ln'] = NEOOutput(np.array(self.a_over_Ln),
                                                   'a/Ln=a(-dln(n0)/dr)',
                                        'normalized equilibrium-scale density')
        self.control['a_over_LT'] = NEOOutput(np.array(self.a_over_LT),
                                                   'a/LT=a(-dln(T0)/dr)',
                                    'normalized equilibrium-scale temperature')
        self.control['inv_tau_self'] = NEOOutput(np.array(self.inv_tau_self),
                                                 'tau_self^-1/(vnorm/a)',
                                         'normalized self-collision frequency')
        self.HH_theory['GAMMA'] = NEOOutput(np.array(self.HH_GAMMA),
                                            'GAMMA/(nnorm*vnorm)',
              'Hinton-Hazeltine second-order radial particle flux (ambipolar)')
        self.HH_theory['Qi'] = NEOOutput(np.array(self.HH_Qi), 'Qi/(nnorm*vnorm*Tnorm)',
                      'Hinton-Hazeltine second-order radial energy flux (ion)')
        self.HH_theory['Qe'] = NEOOutput(np.array(self.HH_Qe), 'Qe/(nnorm*vnorm*Tnorm)',
                 'Hinton-Hazeltine second-order radial energy flux (electron)')
        self.HH_theory['jboot'] = NEOOutput(np.array(self.HH_jboot),
                                            'jpar*B/(e*nnorm*vnorm*Bunit)',
                              'Hinton-Hazeltine first-order bootstrap current')
        self.HH_theory['ki'] = NEOOutput(np.array(self.HH_ki), 'ki',
           'Hinton-Hazeltine first-order dimensionless flow coefficient (ion)')
        self.HH_theory['uipar'] = NEOOutput(np.array(self.HH_uipar),
                                            'uipar*B/(vnorm*Bunit)',
                            'Hinton-Hazeltine first-order parallel flow (ion)')
        self.HH_theory['vipol0'] = NEOOutput(np.array(self.HH_vipol0),
                                             '(vitheta at theta=0)/vnorm',
   'Hinton-Hazeltine first-order poloidal flow at the outboard midplane (ion)')
        self.CH_theory['Qi'] = NEOOutput(np.array(self.CH_Qi), 'Qi/(nnorm*vnorm*Tnorm)',
                          'Chang-Hinton second-order radial energy flux (ion)')
        self.TG_theory['Qi'] = NEOOutput(np.array(self.TG_Qi), 'Qi/(nnorm*vnorm*Tnorm)',
                               'Taguchi second-order radial energy flux (ion)')
        self.S_theory['jboot'] = NEOOutput(np.array(self.S_jboot),
                                           'jpar*B/(e*nnorm*vnorm*Bunit)',
                                        'Sauter first-order bootstrap current')
        self.S_theory['ki'] = NEOOutput(np.array(self.S_ki), 'ki',
                     'Sauter first-order dimensionless flow coefficient (ion)')
        self.S_theory['uipar'] = NEOOutput(np.array(self.S_uipar),
                                           'uipar*B/(vnorm*Bunit)',
                                      'Sauter first-order parallel flow (ion)')
        self.S_theory['vipol0'] = NEOOutput(np.array(self.S_vipol0),
                                            '(vitheta at theta=0)/vnorm',
             'Sauter first-order poloidal flow at the outboard midplane (ion)')
        self.HR_theory['PHI'] = NEOOutput(np.array(self.HR_PHI), '(e*PHI1/Tnorm)^2',
                       'Hinton-Rosenbluth first-order electrostatic potential')
        self.HS_theory['GAMMA'] = NEOOutput(np.array(self.HS_GAMMA),
                                            'GAMMA/(nnorm*vnorm)',
                           'Hirshman-Sigmar second-order radial particle flux')
        self.HS_theory['Q'] = NEOOutput(np.array(self.HS_Q), 'Q/(nnorm*vnorm*Tnorm)',
                             'Hirshman-Sigmar second-order radial energy flux')
        self.transport['GAMMA_gv'] = NEOOutput(np.array(self.SIM_GAMMA_gv),
                                               'GAMMA_gv/(nnorm*vnorm)',
                              'Gyro-viscous second-order radial particle flux')
        self.transport['Q_gv'] = NEOOutput(np.array(self.SIM_Q_gv),
                                           'Q_gv/(nnorm*vnorm*Tnorm)',
                                'Gyro-viscous second-order radial energy flux')
        self.transport['PI_gv'] = NEOOutput(np.array(self.SIM_PI_gv),
                                            'PI_gv/(nnorm*a*Tnorm)',
                              'Gyro-viscous second-order radial momentum flux')

    def read_equil(self):
        """Read out.neo.equil."""

        import sys
        import numpy as np
 
        try:
            equil = np.loadtxt(self.dirname + '/out.neo.equil')
        except:
            print "ERROR (NEOData): Fatal error!  Missing out.neo.equil."
            sys.exit()
            
        for i in range(int(self.n_radial[-1])):
            self.r.append(equil[i, 0])
            self.dPHI0dr.append(equil[i, 1])
            self.q.append(equil[i, 2])
            self.rho_star.append(equil[i, 3])
            self.rmaj_over_a.append(equil[i, 4])
            self.omega_rot.append(equil[i, 5])
            self.omega_rot_deriv.append(equil[i, 6])
            self.nnorm.append(equil[i, 7])
            self.Tnorm.append(equil[i, 8])
            self.vnorm_over_a.append(equil[i, 9])
            temp1 = []
            temp2 = []
            temp3 = []
            temp4 = []
            temp5 = []
            for a in range(int(self.n_species[-1])):
                temp1.append(equil[i, 10+5*a])
                temp2.append(equil[i, 11+5*a])
                temp3.append(equil[i, 12+5*a])
                temp4.append(equil[i, 13+5*a])
                temp5.append(equil[i, 14+5*a])
            self.n0_over_nnorm.append(temp1)
            self.T0_over_Tnorm.append(temp2)
            self.a_over_Ln.append(temp3)
            self.a_over_LT.append(temp4)
            self.inv_tau_self.append(temp5)


    def read_theory(self):
        """Reads out.neo.theory."""

        theory = self.read_file('theory')
        for i in range(int(self.n_radial[-1])):
            self.HH_GAMMA.append(theory[i, 1])
            self.HH_Qi.append(theory[i, 2])
            self.HH_Qe.append(theory[i, 3])
            self.HH_jboot.append(theory[i, 4])
            self.HH_ki.append(theory[i, 5])
            self.HH_uipar.append(theory[i, 6])
            self.HH_vipol0.append(theory[i, 7])
            self.CH_Qi.append(theory[i, 8])
            self.TG_Qi.append(theory[i, 9])
            self.S_jboot.append(theory[i, 10])
            self.S_ki.append(theory[i, 11])
            self.S_uipar.append(theory[i, 12])
            self.S_vipol0.append(theory[i, 13])
            self.HR_PHI.append(theory[i, 14])
            temp1 = []
            temp2 = []
            for a in range(int(self.n_species[-1])):
                temp1.append(theory[i, 15+2*a])
                temp2.append(theory[i, 16+2*a])
            self.HS_GAMMA.append(temp1)
            self.HS_Q.append(temp2)

    def read_transport(self):
        """Reads out.neo.transport."""

        transport = self.read_file('transport')
        for i in range(int(self.n_radial[-1])):
            self.SIM_PHI.append(transport[i, 1])
            self.SIM_jboot.append(transport[i, 2])
            self.SIM_0vtor0.append(transport[i, 3])
            self.SIM_0upar.append(transport[i, 4])
            temp1 = []
            temp2 = []
            temp3 = []
            temp4 = []
            temp5 = []
            temp6 = []
            temp7 = []
            temp8 = []
            for a in range(int(self.n_species[-1])):
                temp1.append(transport[i, 5+8*a])
                temp2.append(transport[i, 6+8*a])
                temp3.append(transport[i, 7+8*a])
                temp4.append(transport[i, 8+8*a])
                temp5.append(transport[i, 9+8*a])
                temp6.append(transport[i, 10+8*a])
                temp7.append(transport[i, 11+8*a])
                temp8.append(transport[i, 12+8*a])
            self.SIM_GAMMA.append(temp1)
            self.SIM_Q.append(temp2)
            self.SIM_PI.append(temp3)
            self.SIM_1upar.append(temp4)
            self.SIM_k.append(temp5)
            self.SIM_K.append(temp6)
            self.SIM_vpol0.append(temp7)
            self.SIM_1vtor0.append(temp8)

    
    def read_transport_gv(self):
        """Reads out.neo.transport_gv."""
        
        transport_gv = self.read_file('transport_gv')
        for i in range(int(self.n_radial[-1])):
            temp1 = []
            temp2 = []
            temp3 = []
            for a in range(int(self.n_species[-1])):
                temp1.append(transport_gv[i, 1+3*a])
                temp2.append(transport_gv[i, 2+3*a])
                temp3.append(transport_gv[i, 3+3*a])
            self.SIM_GAMMA_gv.append(temp1)
            self.SIM_Q_gv.append(temp2)
            self.SIM_PI_gv.append(temp3)

    #-----------------------------------------------#
    # Get data back

    def print_vars(self):
        """Prints all available simulated variables."""

        for key in sorted(self.transport.keys()):
            print key

    #-----------------------------------------------#
    # Plotting routines

    def plot(self, var, n1=2, n2=2, plotcounter=0, fignum=0, legend=True, verbose=False, cols='bgkcmyrw', styles=['-','--','-.',':']):
        """Plots var as a line plot with data from different directories
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

        if plotcounter == 0:
            plotcounter = self.plotcounter
        if fignum == 0:
            fignum = self.fignum

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
        if plotcounter > (n1 * n2):
            plotcounter = 1
            fignum = fignum + 1
        fig = plt.figure(fignum)
        ax = fig.add_subplot(n1, n2, plotcounter)
        ax.set_xlabel('r/a')
        ax.set_ylabel(self.transport.get(var).units)
        #verbose sets title to be more descriptive
        if verbose:
            ax.set_title(self.transport.get(var).descriptor+' vs. r/a')
        else:
            ax.set_title(var + ' vs. r/a')
        tempr = self.control['r'].data
        for var in varlist:

            #If the requested variable is available
            if self.transport.get(var) != None:
                transport = self.transport.get(var).data
                #Check to see if it has one or two dimensions
                try:
                    transport[0][0]
                    #If we reach this point, that means it has two dimensions,
                    #so we have to sort that out properly
                    for y in range(len(transport[0])):
                        plot = []
                        for x in range(len(transport)):
                            plot.append(transport[x][y])
                        plot = np.array(plot)
                        #Create plot of sorted data
                        ax.plot(tempr,plot,c=cols[0],ls=styles[y],label='Sim spe '+str(y))
                except IndexError:
                    #If we reach this point, that means that there is only one
                    #dimension.  All we have to do is create the array, sort it,
                    #and plot it.
                    transport = np.array(transport)
                    ax.plot(tempr, transport, c=cols[0], label='Sim')

            if self.HH_theory.get(var) != None:
                #HH_theory always has only one dimension
                HH = self.HH_theory.get(var).data
                ax.plot(tempr, HH, c=cols[1], label='HH ' + var)

            if self.CH_theory.get(var) != None:
                #So does CH_theory
                CH = self.CH_theory.get(var).data
                ax.plot(tempr, CH, c=cols[2], label='CH ' + var)

            if self.TG_theory.get(var) != None:
                #And TG_theory
                TG = self.TG_theory.get(var).data
                ax.plot(tempr, TG, c=cols[3], label='TG ' + var)

            if self.S_theory.get(var) != None:
                #etc
                S = self.S_theory.get(var).data
                ax.plot(tempr, S, c=cols[4], label='S ' + var)

            if self.HR_theory.get(var) != None:
                HR = self.HR_theory.get(var).data
                ax.plot(tempr, HR, c=cols[5], label='HR ' + var)

            if self.HS_theory.get(var) != None:
                #HS_theory always has two dimensions, so we have to sort that
                #out the same way we sorted out transport
                HS = self.HS_theory.get(var).data
                for y in range(len(HS[0])):
                    plot = []
                    for x in range(len(HS)):
                        plot.append(HS[x][y])
                    plot = np.array(plot)
                    ax.plot(tempr,plot,c=cols[6],ls=styles[y],label='HS '+var+' spe '+str(y))
        #If the legend has been requested, then we create it
        if legend == True:
            ax.legend(loc=2,bbox_to_anchor=(1,1))
        #Move on to the next plot
        plotcounter = plotcounter + 1
        self.plotcounter = plotcounter
        self.fignum = fignum

    #----------------------------------------------------------#
    # Misc methods

    def split(self, array, axis=None):
        """split splits a 2-D array which may be made up of elements with
        multiple entries into an array with one entry per element."""

        import numpy as np


        t = []
        for s in array:
            print s
            if len(s) > 1:
                for a in s:
                    t.append(np.array([a]))
            else:
                t.append(s)
        return t
