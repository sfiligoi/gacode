import data
import sys
import numpy as np
import matplotlib.pyplot as plt
from gacodefuncs import *
from gyro.data import GYROData

class gyrodata_plot(data.GYROData):

    def plot_freq(self,w=0.5,fig=None):
        '''
        Plot gamma and omega vs time

        ARGUMENTS:
         w: fractional width of time window
        '''

        if fig is None:
            fig = plt.figure(figsize=(12,6))

        t = self.t['(c_s/a)t']

        # Read freq data
        self.read_freq()

        # Determine tmin
        for i in range(len(t)):
            if t[i] < (1.0-w)*t[len(t)-1]:
                imin = i

        color = ['k','m','b','c']
        tor_n = self.profile['n0'] + \
                self.profile['d_n']*np.arange(0,self.profile['n_n'])
        #======================================
        ax = fig.add_subplot(121)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(r'$(c_s/a) t$')
        ax.set_ylabel(r'$(a/c_s)\gamma$',color='k')
        #=====================================

        # Gamma
        for i in range(self.profile['n_n']):
            ax.plot(t[imin:],self.freq['(a/c_s)gamma'][i,imin:],color=color[i],
                    label='gamma_n%d'%tor_n[i])

        #======================================
        ax = fig.add_subplot(122)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(r'$(c_s/a) t$')
        ax.set_ylabel(r'$(a/c_s)\omega$',color='k')
        #=====================================

        # Omega
        for i in range(self.profile['n_n']):
            ax.plot(t[imin:],self.freq['(a/c_s)w'][i,imin:],color=color[i],
                    label='omega_n%d'%tor_n[i])

    def plot_balloon(self,index='phi',tmax=-1,fig=None):
        '''
        Plot the ballooning eigenmode structure
       
        ARGUMENTS:
         index: one of 'phi','a','aperp' or an integer
         tmax : max value of (extended angle)/pi
        '''

        if fig is None:
            fig = plt.figure(figsize=(10,6))

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
            ax.set_xlim([-(n_p-1),n_p-1])
        else:
            ax.set_xlim([-tmax,tmax])

        ax.legend()
        return key

    def plot_zf(self,w=0.5,fig=None):
        '''
        Plot the zonal (n=0) potential versus time.

        ARGUMENTS:
         w: fractional time-average window
        '''

        if fig is None:
            fig = plt.figure(figsize=(10,6))
        fig.subplots_adjust(left=0.1,right=0.95,top=0.95,bottom=0.12)

        # Read freq data
        self.read_moment_u()

        ntheta = self.profile['n_theta_plot']
        nx = self.profile['n_x']

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
        print 'INFO: (plot_zf) (nx,ntheta)=',nx,ntheta

        ave_vec = ave*np.ones(len(t))
        #----------------------------------------------------

        ax = fig.add_subplot(111)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(r'$(c_s/a)\, t$')
        ax.set_ylabel(r'$\Phi/\Phi_0$')

        ax.plot(t,y,color='k')
        ax.plot(t[imin:],ave_vec[imin:],color='b',
                label=r'$\mathrm{Average}$',linewidth=1)

        theory = 1.0/(1.0+1.6*2.0**2/np.sqrt(0.5/3.0))
        ax.plot([0,max(t)],[theory,theory],color='grey',
                label=r'$\mathrm{RH \; theory}$',alpha=0.3,linewidth=4)

        ax.legend()
        return ax

    def plot_phi_n0(self,lx=10,ly=6,ymax='auto',span1=-1.0,span2=-1.0,fig=None):
        '''
        Plot the n=0 AND n>0 potentials versus time.

        ARGUMENTS:        
         lx   : width of figure 
         ly   : height of figure 
         ymax : max vertical plot range
         span1: left end of axvspan
         span2: right end of avxspan
        '''
        if fig is None:
            fig = plt.figure(figsize=(lx,ly))
        fig.subplots_adjust(left=0.1,right=0.96,top=0.93,bottom=0.13)

        t = self.t['(c_s/a)t']

        # Read data in gbflux_i and make gbflux
        self.read_field_rms()

        #======================================
        ax = fig.add_subplot(111)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(r'$(c_s/a) t$')
        ax.set_ylabel(r'$\langle e \phi/T_e \rangle/\rho_\star $',color='k')
        #=====================================

        # rho_* for normalization
        x=np.average(self.profile['rho_s'])

        ax.plot(t,self.field_rms[0,:]/x,color='k',label=r'$n=0$')
        ax.plot(t,self.field_rms[1,:]/x,color='purple',label=r'$n>0$')

        if span1 > 0.0:
           ax.axvspan(span1,span2,facecolor='g',alpha=0.1)

        ax.set_xlim([0,t[-1]])

        if ymax != 'auto':
           ax.set_ylim([0,float(ymax)])

        ax.legend()
