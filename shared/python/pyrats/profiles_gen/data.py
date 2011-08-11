"""This data.py contains data classes for input.profiles.

    Contents:
        profiles_genData

"""

class profiles_genData:
    """A class of profiles_gen data.

     Data:

     data
     n_exp
     hlen
     fignum
     plotcounter
     ar
     br
     az
     bz

     Example Usage:
     >>> import matplotlib.pyplot as plt
     >>> from pyrats.profiles_gen.data import profiles_genData
     >>> prof1 = profiles_genData('$GACODE_ROOT/tgyro/tools/input/treg01')
     >>> prof1.plot('ne')
     >>> plt.show()

"""


    #Methods
    def __init__(self, directory = '.'):
        """Constructor reads in data from directory and creates new object."""
        
        self.set_directory(directory)
        self.init_data()
        self.store_data()

    def init_data(self):
        """Initialize object data."""

        self.data = {}
        self.n_exp = 0
        self.hlen = 0
        self.fignum = 1
        self.plotcounter = 1
        self.ar = []
        self.br = []
        self.az = []
        self.bz = []

    def set_directory(self, directory):
        """Set the directory which contains input.profiles."""

        from os.path import expanduser, expandvars
        self.directory_name = expanduser(expandvars(directory))

    def read_data(self):
        """Read in object data from input.profiles."""

        import numpy as np

        elements = {}
        temp = []
        raw_data = open(self.directory_name + '/input.profiles', 'r').readlines()
        #Determines length of header, and reads in data
        while raw_data[self.hlen].strip()[0].isdigit() == False:
            self.hlen = self.hlen + 1
            if raw_data[self.hlen].strip()[0:6] == 'N_EXP=':
                self.n_exp = int(raw_data[self.hlen].strip()[6:])
        for line in range(self.hlen, len(raw_data)):
            if raw_data[line].strip()[0].isdigit() or raw_data[line].strip()[0] == '-':
                temp.append(raw_data[line].split())
        self.hlen = self.hlen - 1
        data = np.array(temp)

        #Reads in variable names
        keywords = []
        for count in range(5):
            #This line separates columns when they run together
            keywords = raw_data[(self.n_exp + 2)*count + self.hlen].replace('(kW/eV)', '(kW/eV) ').split()
            column = 0
            for key in keywords:
                elements[key] = data[0:self.n_exp, column]
                column = column + 1
            data = data[self.n_exp:, :]
        #We have to separate the keys into two chunks because the spacing
        #changes after the 5th group of variables, and so we need a different
        #method for reading them.
        for count in range(5, 8):
            keys = raw_data[(self.n_exp + 2) * count + self.hlen].split('  ')
            c = keys.count('')
            for n in range(c):
                keys.remove('')
            keywords = []
            for x in range(5):
                keywords.append(keys[x])
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
        """Reads data and renames it appropriately.

        store_data is necessary because the names of the different parameters are
        not uniformly formatted.  store_data cleans them up by inserting spaces
        where necessary, and by deleting #-signs when necessary."""

        self.data = self.read_data()
        for k in self.data.keys():
            flag = False
            temp = list(k)
            for i in range(len(temp)):
                if temp[i] == '#':
                    del temp[i]
                    flag = True
                    break
            for i in range(len(temp)):
                if temp[i] == '(':
                    if temp[i-1] != ' ':
                        temp.insert(i, ' ')
                        flag = True
                        break
            if flag:
                self.data[''.join(temp)] = self.data[k]
                self.data.pop(k)

    def compute_mtypeeq(self, r):

        """Uses input data to compute the Miller-type equilibrium at radius r."""
        
        import math

        R = []
        Z = []
        tot = []
        theta = 0
        dtheta = 0.01
        x = self.match(r, self.data['rho (-)'])
        while theta < 2 * math.pi:
            a = float(self.data['rmaj (m)'][x]) + float(self.data['rmin (m)'][x]) * math.cos(theta + math.asin(float(self.data['delta (-)'][x])) * math.sin(theta))
            R.append(a)
            b = float(self.data['zmag (m)'][x]) + float(self.data['kappa (-)'][x]) * float(self.data['rmin (m)'][x]) * math.sin(theta + float(self.data['zeta (-)'][x]) * math.sin(2 * theta))
            Z.append(b)
            theta = theta + dtheta
        tot.append(R)
        tot.append(Z)
        tot.append(self.data['rmaj (m)'][x])
        tot.append(self.data['zmag (m)'][x])
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
        x = self.match(r, self.data['rho (-)'])
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
        tot.append(self.data['rmaj (m)'][x])
        tot.append(self.data['zmag (m)'][x])
        return tot

    #-------------------------------------------- #
    # Get data back
    def get(self, var):
        """Return requested variable."""

        return self.data[var]

    #-------------------------------------------- #
    # Plotting functions
    def plot(self, var, n1=2, n2=2, plotcounter=0, fignum=0):
        """Plots requested data using matplotlib.

        var is the requested variable to be plotted,
        n1 is the horizontal number of plots in one window.
        n2 is the vertical number of plots in one window.

        Examples: self.plot(rmaj, 2, 2)"""

        import matplotlib.pyplot as plt
        import matplotlib as mpl

        mpl.rcParams['figure.subplot.wspace'] = .3
        mpl.rcParams['figure.subplot.hspace'] = .4

        if plotcounter == 0:
            plotcounter = self.plotcounter
        if fignum == 0:
            fignum = self.fignum

        if plotcounter > (n1 * n2):
            plotcounter = 1
            fignum = fignum + 1
        fig = plt.figure(fignum)
        ax = fig.add_subplot(n2, n1, plotcounter)
        for k in self.data.iterkeys():
            if (var + ' ') == k[:len(var) + 1]:
                toplot = k
                ylab = k
                ax.set_xlabel(k)
        ax.set_xlabel(u'\u03c1 (-)')
        if var == 'kappa':
            ylab = u'\u03ba (-)'
        if var == 'delta':
            ylab = u'\u03b4 (-)'
        if var == 'zeta':
            ylab = u'\u03b6 (-)'
        if var == 'omega0':
            ylab = u'\u03c90 (1/s)'
        ax.set_ylabel(ylab)
        ax.set_title(ylab.split()[0] + u' vs. \u03c1')
        ax.plot(self.data['rho (-)'], self.data[toplot])
        plotcounter = plotcounter + 1
        self.plotcounter = plotcounter
        self.fignum = fignum


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
        fig = plt.figure(self.fignum)
        self.fignum = self.fignum + 1
        ax = fig.add_subplot(111, aspect='equal')
        ax.set_ylabel('Z (m)')
        ax.set_xlabel('R (m)')
        ax.axhline(c='k', ls='--')
        #Increments through from min to max and creates plots at each
        #radius.
        inc = self.match(inner, self.data['rho (-)'])
        step = (self.match(outer, self.data['rho (-)']) - inc)/n
        count = 0
        while inc < self.match(outer, self.data['rho (-)']):
            mteq.append(self.compute_mtypeeq(float(self.data['rho (-)'][inc])))
            r = mteq[0][0]
            z = mteq[0][1]
            rmaj = float(mteq[0][2])
            zmag = float(mteq.pop()[3])
            ax.plot(r, z, 'b', label='Miller surface at '+str(self.data['rho (-)'][inc]))
            count = count + math.modf(step)[0]
            if count > 1:
                count = math.modf(count)[0]
                inc = inc + math.modf(step)[1] + 1
            else:
                inc = inc + math.modf(step)[1]
        #Sets limits, just as above.
        if n == 1:
            ax.set_title(u'Flux Surface at \u03c1 = ' + str(self.data['rho (-)'][self.match(inner, self.data['rho (-)'])]))
            ax.axhline(y=self.get('zmag (m)')[self.match(inner,self.get('rho (-)'))], c='m')
            ax.axvline(x=self.get('rmaj (m)')[self.match(inner, self.get('rho (-)'))], c='m')
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
        fig = plt.figure(self.fignum)
        self.fignum = self.fignum + 1
        ax = fig.add_subplot(111, aspect='equal')
        ax.set_ylabel('Z (m)')
        ax.set_xlabel('R (m)')
        ax.axhline(c='k', ls='--')
        inc = self.match(inner, self.data['rho (-)'])
        step = (self.match(outer, self.data['rho (-)']) - inc)/n
        count = 0
        while inc < self.match(outer, self.data['rho (-)']):
            feq.append(self.compute_fouriereq(float(self.data['rho (-)'][inc])))
            r = feq[0][0]
            z = feq[0][1]
            rmaj = float(feq[0][2])
            zmag = float(feq.pop()[3])
            lines.append(ax.plot(r, z, 'r', label='Fourier surface at '+str(self.data['rho (-)'][inc])))
            count = count + math.modf(step)[0]
            if count > 1:
                count = math.modf(count)[0]
                inc = inc + math.modf(step)[1] + 1
            else:
                inc = inc + math.modf(step)[1]
        if n == 1:
            ax.set_title(u'Flux Surface at \u03c1 = ' + str(self.data['rho (-)'][self.match(inner, self.data['rho (-)'])]))
            ax.axhline(y=self.get('zmag (m)')[self.match(inner,self.get('rho (-)'))], c='m')
            ax.axvline(x=self.get('rmaj (m)')[self.match(inner, self.get('rho (-)'))], c='m')
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
        fig = plt.figure(self.fignum)
        self.fignum = self.fignum + 1
        ax = fig.add_subplot(111, aspect='equal')
        ax.set_ylabel('Z (m)')
        ax.set_xlabel('R (m)')
        ax.axhline(c='k', ls='--')
        inc = self.match(inner, self.data['rho (-)'])
        step = (self.match(outer, self.data['rho (-)']) - inc)/n
        count = 0
        while inc < self.match(outer, self.data['rho (-)']):
            feq.append(self.compute_fouriereq(float(self.data['rho (-)'][inc])))
            fr = feq[0][0]
            fz = feq.pop()[1]
            ax.plot(fr, fz, 'r', label='Fourier surface at '+str(self.data['rho (-)'][inc]))
            mteq.append(self.compute_mtypeeq(float(self.data['rho (-)'][inc])))
            mr = mteq[0][0]
            mz = mteq[0][1]
            rmaj = float(mteq[0][2])
            zmag = float(mteq.pop()[3])
            ax.plot(mr, mz, 'b', label='Miller surface at '+str(self.data['rho (-)'][inc]))
            count = count + math.modf(step)[0]
            if count > 1:
                count = math.modf(count)[0]
                inc = inc + math.modf(step)[1] + 1
            else:
                inc = inc + math.modf(step)[1]
        if n == 1:
            ax.set_title(u'Flux Surface at \u03c1 = ' + str(self.data['rho (-)'][self.match(inner, self.data['rho (-)'])]))
            ax.axhline(y=self.get('zmag (m)')[self.match(inner,self.get('rho (-)'))], c='m')
            ax.axvline(x=self.get('rmaj (m)')[self.match(inner, self.get('rho (-)'))], c='m')
            ax.legend( ('Midplane', 'Miller-type', 'Fourier-type', 'Flux surface center'), loc=2,bbox_to_anchor=(1,1))
        elif verbose:
            ax.set_title(str(int(n)) + u' Flux Surfaces between \u03c1 = ' + str(inner) + u' and \u03c1 = ' + str(outer) + '.')
            ax.legend(loc=2,bbox_to_anchor=(1,1))
        else:
            ax.set_title(str(int(n)) + u' Flux Surfaces between \u03c1 = ' + str(inner) + u' and \u03c1 = ' + str(outer) + '.')
            ax.legend( ('Midplane', 'Miller-type', 'Fourier-type'), loc=2,bbox_to_anchor=(1,1))
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
