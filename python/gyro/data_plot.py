import data
import sys
import numpy as np
from gacodeplotdefs import *
from gyro.data import GYROData

#---------------------------------------------------------------------------#
class gyrodata_plot(data.GYROData):
    def plot_freq(self, window=0.5, fig=None):
        """Plot frequency vs time"""
        """window = avg. window fraction"""

        GFONTSIZE=16
        if (fig is None):
            fig = plt.figure()

        t    = self.t['(c_s/a)t']

        # Read freq data
        self.read_freq()

        # Determine tmin
        for i in range(len(t)):
            if t[i] < (1.0-window)*t[len(t)-1]:
                imin = i

        color = ['k','m','b','c']
        tor_n = self.profile['n0'] + \
                self.profile['d_n']*np.arange(0,self.profile['n_n'])
        #======================================
        ax = fig.add_subplot(121)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(r'$(c_s/a) t$',fontsize=GFONTSIZE)
        ax.set_ylabel(r'$(a/c_s)\gamma$',color='k',fontsize=GFONTSIZE)
        #=====================================

        #Gamma
        for i in range(self.profile['n_n']):
            ax.plot(t[imin:],self.freq['(a/c_s)gamma'][i,imin:],color=color[i],
                    label='gamma_n%d'%tor_n[i])

        #======================================
        ax = fig.add_subplot(122)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(r'$(c_s/a) t$',fontsize=GFONTSIZE)
        ax.set_ylabel(r'$(a/c_s)\omega$',color='k',fontsize=GFONTSIZE)
        #=====================================

        #Omega
        for i in range(self.profile['n_n']):
            ax.plot(t[imin:],self.freq['(a/c_s)w'][i,imin:],color=color[i],
                    label='omega_n%d'%tor_n[i])

    def plot_balloon(self, index='phi',tmax=-1, fig=None):
        '''
        Plot the ballooning mode structure

        index - One of 'phi','a','aperp' or an integer
        '''
        self.read_balloon()

        if index == 'phi':
            key = 'balloon_phi'
            index = self.balloon.keys().index(key)
        elif index == 'a':
            key = 'balloon_a'
            index = self.balloon.keys().index(key)
        elif index == 'aperp':
            key = 'balloon_aperp'
            index = self.balloon.keys().index(key)
        else:
            index = int(index)

        key = self.balloon.keys()[index]

        ytitle = '\\begin{verbatim}'+key+'\\end{verbatim}'

        if key == 'balloon_a':
            ytitle = r'$\delta A_\parallel$'
        if key == 'balloon_phi':
            ytitle = r'$\delta \phi$'
        if key == 'balloon_aperp':
            ytitle = r'$\delta B_\parallel$'
        if key == 'balloon_epar':
            ytitle = r'$\delta E_\parallel$'

        #======================================
        if fig is None:
            fig = plt.figure(figsize=(10,6))
        ax = fig.add_subplot(111)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(r'$\theta_*/\pi$')
        ax.set_ylabel(ytitle)
        #=====================================

        n_p   = self.profile['n_x']/self.profile['box_multiplier']
        n_ang = self.profile['n_theta_plot']*n_p

        x = -(1.0+n_p)+2.0*n_p*np.arange(n_ang)/float(n_ang)

        ax.plot(x,np.real(self.balloon[key][:,0,-1]),color='k',label='Re')
        ax.plot(x,np.imag(self.balloon[key][:,0,-1]),color='m',label='Im')

        if tmax < 0.0:
            ax.set_xlim([1-n_p,-1+n_p])
        else:
            ax.set_xlim([-tmax,tmax])

        ax.legend()
        return key

    def plot_zf(self,w=.5,fig=None):
        '''
        plot the zonal flow?

        w - how much of normalized time to use
        '''

        # Read freq data
        self.read_moment_u()

        ntheta = self.profile['n_theta_plot']
        nx = self.profile['n_x']

        print nx,ntheta
        y = np.real(self.moment_u[ntheta/3,nx/3,0,0,:])

        y = y[:]/y[0]
        t = self.t['(c_s/a)t']

        #----------------------------------------------------
        # Average calculations

        # Determine tmin
        for i in range(len(t)):
            if t[i] < (1.0-w)*t[len(t)-1]:
                imin = i

        ave = average(y[:],t,w)
        print 'INFO: (plot_zf) Integral time-average = %.6f' % ave

        ave_vec = ave*np.ones(len(t))
        #----------------------------------------------------

        #======================================
        if fig is None:
            fig = plt.figure(figsize=(10,6))
        fig.subplots_adjust(left=0.1,right=0.95,top=0.95,bottom=0.12)

        ax = fig.add_subplot(111)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(r'$(c_s/a)\, t$',fontsize=GFONTSIZE)
        ax.set_ylabel(r'$\Phi/\Phi_0$',fontsize=GFONTSIZE)

        ax.plot(t,y,color='k')
        ax.plot(t[imin:],ave_vec[imin:],color='b',label=r'$\mathrm{Average}$',linewidth=1)

        theory = 1.0/(1.0+1.6*2.0**2/np.sqrt(0.5/3.0))
        ax.plot([0,max(t)],[theory,theory],color='grey',label=r'$\mathrm{RH \; theory}$',alpha=0.3,linewidth=4)

        ax.legend()
