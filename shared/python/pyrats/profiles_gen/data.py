class profiles_genData:
    """profiles_gen output data (i.e., input.profiles) class.

     Data:

     data
     n_exp
     hlen
     ar
     br
     az
     bz
     max_rho

     Example Usage:
     >>> import matplotlib.pyplot as plt
     >>> from pyrats.profiles_gen.data import profiles_genData
     >>> prof = profiles_genData('input.profiles')

     """

    #---------------------------------------------------------------------------#
    # Methods

    def __init__(self,infile):
        """
        Constructor reads in data from directory and creates new object.
        """

        import string
        import numpy as np

        self.infile=infile

        # Initialize data
        self.data = {}
        self.n_exp = 0
        self.hlen = 0
        self.ar = []
        self.br = []
        self.az = []
        self.bz = []

        # Read data
        row  = 0
        for line in open(infile,'r').readlines():
            row = row+1
            if line[0:5] == 'N_EXP':
                self.n_exp = int(string.splitfields(line,'=')[1])
            if line[0:8] == 'ARHO_EXP':
                self.max_rho = float(line.split('=')[1])
            if line[0:7] == '#rho(-)':
                break

        data = np.loadtxt(infile,skiprows=row)

        n = self.n_exp

        self.data['rho']   = data[0:n,0]
        self.data['rmin']  = data[0:n,1]
        self.data['rmaj']  = data[0:n,2]
        self.data['q']     = data[0:n,3]
        self.data['kappa'] = data[0:n,4]

        self.data['delta']  = data[n:2*n,0]
        self.data['Te']     = data[n:2*n,1]
        self.data['ne']     = data[n:2*n,2]
        self.data['z_eff']  = data[n:2*n,3]
        self.data['omega0'] = data[n:2*n,4]

        self.data['flow_mom'] = data[2*n:3*n,0]
        self.data['pow_e']    = data[2*n:3*n,1]
        self.data['pow_i']    = data[2*n:3*n,2]
        self.data['pow_ei']   = data[2*n:3*n,3]
        self.data['zeta']     = data[2*n:3*n,4]

        self.data['flow_beam'] = data[3*n:4*n,0]
        self.data['flow_wall'] = data[3*n:4*n,1]
        self.data['zmag']      = data[3*n:4*n,2]
        self.data['ptot']      = data[3*n:4*n,3]
        self.data['polflux']   = data[3*n:4*n,4]

        self.data['ni_1'] = data[4*n:5*n,0]
        self.data['ni_2'] = data[4*n:5*n,1]
        self.data['ni_3'] = data[4*n:5*n,2]
        self.data['ni_4'] = data[4*n:5*n,3]
        self.data['ni_5'] = data[4*n:5*n,4]

        self.data['Ti_1'] = data[5*n:6*n,0]
        self.data['Ti_2'] = data[5*n:6*n,1]
        self.data['Ti_3'] = data[5*n:6*n,2]
        self.data['Ti_4'] = data[5*n:6*n,3]
        self.data['Ti_5'] = data[5*n:6*n,4]

        self.data['vtor_1'] = data[6*n:7*n,0]
        self.data['vtor_2'] = data[6*n:7*n,1]
        self.data['vtor_3'] = data[6*n:7*n,2]
        self.data['vtor_4'] = data[6*n:7*n,3]
        self.data['vtor_5'] = data[6*n:7*n,4]

        self.data['vpol_1'] = data[7*n:8*n,0]
        self.data['vpol_2'] = data[7*n:8*n,1]
        self.data['vpol_3'] = data[7*n:8*n,2]
        self.data['vpol_4'] = data[7*n:8*n,3]
        self.data['vpol_5'] = data[7*n:8*n,4]

        data = np.loadtxt(infile+'.extra')
 
        try:
            data = np.loadtxt(infile+'.extra',comments='#')
            x = data.reshape((n,35), order='F')

            self.data['bunit']     = x[0:n,0]
            self.data['s']         = x[0:n,1]
            self.data['drmaj']     = x[0:n,2]
            self.data['dzmag']     = x[0:n,3]
            self.data['sdelta']    = x[0:n,4]
            self.data['skappa']    = x[0:n,5]
            self.data['szeta']     = x[0:n,6]
            self.data['dlnnedr']   = x[0:n,7]
            self.data['dlntedr']   = x[0:n,8]
            self.data['dlnnidr_1'] = x[0:n,9]
            self.data['dlnnidr_2'] = x[0:n,10]
            self.data['dlnnidr_3'] = x[0:n,11]
            self.data['dlnnidr_4'] = x[0:n,12]
            self.data['dlntidr_5'] = x[0:n,13]
            self.data['dlntidr_1'] = x[0:n,14]
            self.data['dlntidr_2'] = x[0:n,15]
            self.data['dlntidr_3'] = x[0:n,16]
            self.data['dlntidr_4'] = x[0:n,17]
            self.data['dlntidr_5'] = x[0:n,18]
            self.data['dlnptotdr'] = x[0:n,19]
            self.data['drdrho']    = x[0:n,20]
            self.data['w0p']       = x[0:n,21]
            self.data['vol']       = x[0:n,22]
            self.data['volp']      = x[0:n,23]
            self.data['cs']        = x[0:n,24]
            self.data['rhos']      = x[0:n,25]
            self.data['ni_new']    = x[0:n,26]
            self.data['dlnnidr_new'] = x[0:n,27]
            self.data['grad_r0']   = x[0:n,28]
            self.data['ave_grad_r']= x[0:n,29]
            self.data['bp0']       = x[0:n,30]
            self.data['bt0']       = x[0:n,31]
            self.data['gamma_e']   = x[0:n,32]
            self.data['gamma_p']   = x[0:n,33]
            self.data['mach']      = x[0:n,34]
        
        except:
            print infile+'.extra not available.'


    def read_fourier(self):
        """Read in data from input.profiles.geo."""

        temp = []
        raw_data = open(self.infile+'.geo', 'r').readlines()
        for line in raw_data:
            temp.append(line.split()[0])
        count = int(temp[0]) + 1

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


    def compute_mtypeeq(self, r):

        """Uses input data to compute the Miller-type equilibrium at radius r."""
        
        import math

        R = []
        Z = []
        tot = []
        theta = 0
        dtheta = 0.01
        x = self.match(r, self.data['rho'])
        while theta < 2 * math.pi:
            a = float(self.data['rmaj'][x]) + float(self.data['rmin'][x])*math.cos(theta + math.asin(float(self.data['delta'][x])) * math.sin(theta))
            R.append(a)
            b = float(self.data['zmag'][x]) + float(self.data['kappa'][x])*float(self.data['rmin'][x])*math.sin(theta + float(self.data['zeta'][x])*math.sin(2 * theta))
            Z.append(b)
            theta = theta + dtheta
        tot.append(R)
        tot.append(Z)
        tot.append(self.data['rmaj'][x])
        tot.append(self.data['zmag'][x])
        return tot

    def compute_fouriereq(self, r):

        """Uses input data to compute the general Fourier-series equilibrium at
        radius r."""

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


    def millerplot(self, inner, outer, n, verbose):
        """Creates plots of Miller-type equilibrium flux surfaces.

        This method will plot n flux surfaces between inner and outer.  If
        verbose is set to true, and n > 1, then the legend will display the
        locations of each flux surface.
        """
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import sys
        import math

        mpl.rcParams['font.size'] = 10.0
        mpl.rcParams['figure.subplot.right'] = .7
        mteq = []
        bigcoords = [-100, 100]
        midplane = [0, 0]

        #Produces a matplotlib figure object and creates the labels.
        fig = plt.figure(figsize=(12,7))
        ax = fig.add_subplot(111, aspect='equal')
        ax.set_ylabel('Z (m)')
        ax.set_xlabel('R (m)')
        ax.axhline(c='k', ls='--')
        #Increments through from min to max and creates plots at each
        #radius.
        inc = self.match(inner, self.data['rho'])
        step = (self.match(outer, self.data['rho']) - inc)/n
        count = 0
        while inc < self.match(outer, self.data['rho']):
            mteq.append(self.compute_mtypeeq(float(self.data['rho'][inc])))
            r = mteq[0][0]
            z = mteq[0][1]
            rmaj = float(mteq[0][2])
            zmag = float(mteq.pop()[3])
            ax.plot(r, z, 'b', label='Miller surface at '+str(self.data['rho'][inc]))
            count = count + math.modf(step)[0]
            if count > 1:
                count = math.modf(count)[0]
                inc = inc + math.modf(step)[1] + 1
            else:
                inc = inc + math.modf(step)[1]
        #Sets limits, just as above.
        if n == 1:
            ax.set_title(u'Flux Surface at \u03c1 = ' + str(self.data['rho'][self.match(inner, self.data['rho'])]))
            ax.axhline(y=self.get('zmag')[self.match(inner,self.get('rho'))], c='m')
            ax.axvline(x=self.get('rmaj')[self.match(inner, self.get('rho'))], c='m')
            ax.legend( ('Midplane', 'Miller-type', 'Flux surface center'), loc=2,bbox_to_anchor=(1,1))
        elif verbose:
            ax.set_title(str(int(n)) + u' Flux Surfaces between \u03c1 = ' + str(inner) + u' and \u03c1 = ' + str(outer) + '.')
            ax.legend(loc=2,bbox_to_anchor=(1,1))
        else:
            ax.set_title(str(int(n)) + u' Flux Surfaces between \u03c1 = ' + str(inner) + u' and \u03c1 = ' + str(outer) + '.')
            ax.legend( ('Midplane', 'Miller-type'), loc=2,bbox_to_anchor=(1,1))
        aspect = max((max(z) - min(z)), (max(r) - min(r))) + .5
        ax.set_xlim(rmaj - aspect/2, rmaj + aspect/2)
        ax.set_ylim(zmag - aspect/2, zmag + aspect/2)

    def fourierplot(self, inner, outer, n, verbose):
        """Creates plots of Fourier-type equilibrium flux surfaces.

        This method will plot n flux surfaces between inner and outer.  If
        verbose is set to true, and n > 1, then the legend will display the
        locations of each flux surface.
        """

        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import sys
        import math

        mpl.rcParams['font.size'] = 10.0
        mpl.rcParams['figure.subplot.right'] = .7
        feq = []
        lines = []
        bigcoords = [-100, 100]
        midplane = [0, 0]

        #Produces a matplotlib figure object and creates the labels.
        fig = plt.figure(figsize=(12,7))
        ax = fig.add_subplot(111, aspect='equal')
        ax.set_ylabel('Z (m)')
        ax.set_xlabel('R (m)')
        ax.axhline(c='k', ls='--')
        inc = self.match(inner, self.data['rho'])
        step = (self.match(outer, self.data['rho']) - inc)/n
        count = 0
        while inc < self.match(outer, self.data['rho']):
            feq.append(self.compute_fouriereq(float(self.data['rho'][inc])))
            r = feq[0][0]
            z = feq[0][1]
            rmaj = float(feq[0][2])
            zmag = float(feq.pop()[3])
            lines.append(ax.plot(r, z, 'r', label='Fourier surface at '+str(self.data['rho'][inc])))
            count = count + math.modf(step)[0]
            if count > 1:
                count = math.modf(count)[0]
                inc = inc + math.modf(step)[1] + 1
            else:
                inc = inc + math.modf(step)[1]
        if n == 1:
            ax.set_title(u'Flux Surface at \u03c1 = ' + str(self.data['rho'][self.match(inner, self.data['rho'])]))
            ax.axhline(y=self.get('zmag')[self.match(inner,self.get('rho'))], c='m')
            ax.axvline(x=self.get('rmaj')[self.match(inner, self.get('rho'))], c='m')
            ax.legend( ('Midplane', 'Fourier-type', 'Flux surface center'), loc=2,bbox_to_anchor=(1,1))
        elif verbose:
            ax.set_title(str(int(n)) + u' Flux Surfaces between \u03c1 = ' + str(inner) + u' and \u03c1 = ' + str(outer) + '.')
            ax.legend(loc=2,bbox_to_anchor=(1,1))
        else:
            ax.set_title(str(int(n)) + u' Flux Surfaces between \u03c1 = ' + str(inner) + u' and \u03c1 = ' + str(outer) + '.')
            ax.legend( ('Midplane', 'Fourier-type'), loc=2,bbox_to_anchor=(1,1))
        aspect = max((max(z) - min(z)), (max(r) - min(r))) + .5
        ax.set_xlim(rmaj - aspect/2, rmaj + aspect/2)
        ax.set_ylim(zmag - aspect/2, zmag + aspect/2)

    def compplot(self, inner, outer, n, verbose):
        """Creates plots of both Miller-type and Fourier-type equilibrium flux
        surfaces.

        This method will plot n flux surfaces between inner and outer.  If
        verbose is set to true, and n > 1, then the legend will display the
        locations of each flux surface.
        """

        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import sys
        import math

        mpl.rcParams['font.size'] = 10.0
        mpl.rcParams['figure.subplot.right'] = .7
        mteq = []
        lines = []
        feq = []
        bigcoords = [-100, 100]
        midplane = [0, 0]

        #Produces a matplotlib figure object and creates the labels.
        fig = plt.figure(figsize=(12,7))
        ax = fig.add_subplot(111, aspect='equal')
        ax.set_ylabel('Z (m)')
        ax.set_xlabel('R (m)')
        ax.axhline(c='k', ls='--')
        inc = self.match(inner, self.data['rho'])
        step = (self.match(outer, self.data['rho']) - inc)/n
        count = 0
        while inc < self.match(outer, self.data['rho']):
            feq.append(self.compute_fouriereq(float(self.data['rho'][inc])))
            fr = feq[0][0]
            fz = feq.pop()[1]
            ax.plot(fr, fz, 'r', label='Fourier surface at '+str(self.data['rho'][inc]))
            mteq.append(self.compute_mtypeeq(float(self.data['rho'][inc])))
            mr = mteq[0][0]
            mz = mteq[0][1]
            rmaj = float(mteq[0][2])
            zmag = float(mteq.pop()[3])
            ax.plot(mr, mz, 'b', label='Miller surface at '+str(self.data['rho'][inc]))
            count = count + math.modf(step)[0]
            if count > 1:
                count = math.modf(count)[0]
                inc = inc + math.modf(step)[1] + 1
            else:
                inc = inc + math.modf(step)[1]
        if n == 1:
            ax.set_title(u'Flux Surface at \u03c1 = ' + str(self.data['rho'][self.match(inner, self.data['rho'])]))
            ax.axhline(y=self.get('zmag')[self.match(inner,self.get('rho'))], c='m')
            ax.axvline(x=self.get('rmaj')[self.match(inner, self.get('rho'))], c='m')
            ax.legend( ('Midplane', 'Fourier-type', 'Miller-type', 'Flux surface center'), loc=2,bbox_to_anchor=(1,1))
        elif verbose:
            ax.set_title(str(int(n)) + u' Flux Surfaces between \u03c1 = ' + str(inner) + u' and \u03c1 = ' + str(outer) + '.')
            ax.legend(loc=2,bbox_to_anchor=(1,1))
        else:
            ax.set_title(str(int(n)) + u' Flux Surfaces between \u03c1 = ' + str(inner) + u' and \u03c1 = ' + str(outer) + '.')
            ax.legend( ('Midplane', 'Fourier-type', 'Miller-type'), loc=2,bbox_to_anchor=(1,1))
        aspect = max((max(fz) - min(fz)), (max(fr) - min(fr)), (max(mz) - min(mz)), (max(mr) - min(mr))) + .5
        ax.set_xlim(rmaj - aspect/2, rmaj + aspect/2)
        ax.set_ylim(zmag - aspect/2, zmag + aspect/2)

    #-------------------------------------------- #
    # Misc Functions
    def match(self, val, vec):
        """Return index of closest match to val in a list of values (vec)."""
        import math

        temp = []
        for n in range(len(vec)):
            t = math.fabs(val - float(vec[n])), n
            temp.append(t)
        return sorted(temp)[0][1]
