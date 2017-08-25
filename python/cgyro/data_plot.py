import data
import sys
import numpy as np
import matplotlib.pyplot as plt
from gacodefuncs import *
from cgyro.data import cgyrodata

class cgyrodata_plot(data.cgyrodata):

   def plot_freq(self,w=0.5,fig=None):
      '''
      Plot gamma and omega vs time

      ARGUMENTS:
      w: fractional width of time window
      '''

      if fig is None:
         fig = plt.figure(figsize=(12,6))
         fig.subplots_adjust(left=0.085,right=0.97,top=0.92,bottom=0.12)

      #======================================
      # Omega
      ax = fig.add_subplot(121)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)
      ax.set_ylabel(r'$(a/c_s)\, \omega$')

      for i in range(self.n_n):
         ax.plot(self.t,self.freq[0,i,:])

      ax.set_xlim([0,self.t[-1]])
      #======================================

      #======================================
      # Gamma
      ax = fig.add_subplot(122)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)
      ax.set_ylabel(r'$(a/c_s)\, \gamma$')

      for i in range(self.n_n):
         ax.plot(self.t,self.freq[1,i,:])

      ax.set_xlim([0,self.t[-1]])
      #======================================

      
   def plot_ky_freq(self,w=0.5,fig=None):
      '''
      Plot gamma and omega vs ky

      ARGUMENTS:
      w: fractional width of time window
      '''

      if fig is None:
         fig = plt.figure(figsize=(12,6))
         fig.subplots_adjust(left=0.085,right=0.97,top=0.92,bottom=0.12)

      #======================================
      # Omega
      ax = fig.add_subplot(121)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(r'$k_y \rho_s$')
      ax.set_ylabel(r'$(a/c_s)\, \omega$')

      ax.plot(self.ky,self.freq[0,:,-1],color='blue')
      ax.plot(self.ky,self.freq[0,:,-1],"o",color='k')
      if len(self.ky) > 1:
         ax.set_xlim([0,self.ky[-1]])
      #======================================

      #======================================
      # Gamma
      ax = fig.add_subplot(122)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(r'$k_y \rho_s$')
      ax.set_ylabel(r'$(a/c_s)\, \gamma$')
         
      ax.plot(self.ky,self.freq[1,:,-1],color='red')
      ax.plot(self.ky,self.freq[1,:,-1],"o",color='k')
      if len(self.ky) > 1:
         ax.set_xlim([0,self.ky[-1]])
      #======================================
    
   def plot_ky_phi(self,field=0,ymin='0',ymax='auto',nstr='null',fig=None):
      '''
      Plot gamma and omega vs ky

      ARGUMENTS:
      w: fractional width of time window
      '''

      if fig is None:
         fig = plt.figure(figsize=(12,6))
         fig.subplots_adjust(left=0.08,right=0.96,top=0.92,bottom=0.12)

      self.getbigfield()

      ax = fig.add_subplot(111)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)
      ax.set_ylabel(r'$\left| \delta\phi_n \right|$')
      ax.set_yscale('log')
      ax.set_title(r'$\mathrm{Fluctuation~intensity} \quad k_\theta = nq/r$')

      p2 = np.sum(self.phisq[:,0,:,:],axis=0)/self.rho**2

      if nstr == 'null':
         nvec = range(self.n_n)
      else:
         nvec = str2list(nstr)
    
      for n in nvec:
         num = r'$n='+str(n)+'$'
         if n==0:
            ax.plot(self.t,np.sqrt(p2[n,:]),linewidth=2,label=num)
         else:
            ax.plot(self.t,np.sqrt(p2[n,:]),label=num)

      ax.set_xlim([0,max(self.t)])

      if self.n_n > 16:
         ax.legend(loc=4, ncol=5, prop={'size':12})
      else:
         ax.legend(loc=4, ncol=6, prop={'size':12})

      ymin,ymax=setlimits(ax.get_ylim(),ymin,ymax)

      ax.set_ylim([ymin,ymax])


   def plot_rcorr_phi(self,w=0.5,fig=None):
      '''
      Plot gamma and omega vs ky

      ARGUMENTS:
      w: fractional width of time window
      '''

      import scipy as scipy
      import scipy.signal as signal
      from scipy.optimize import curve_fit

      def absexp(x,tau):
         return np.exp(-np.abs(x)/tau)

      if fig is None:
         fig = plt.figure(figsize=(12,6))
         fig.subplots_adjust(left=0.08,right=0.96,top=0.92,bottom=0.12)

      self.getbigfield()

      t   = self.t
      kx  = self.kx
      ave = np.zeros(self.n_radial)

      imin=iwindow(self.t,w)
    
      dk = kx[1]-kx[0]
      x0 = kx[-1]+dk

      color = ['m','k','b','c']
      xlabel=r'$r / \rho_s$'
      windowtxt = r'$['+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+']$'

      ax = fig.add_subplot(1,1,1)
      ax.set_title(r'$\mathrm{Average~radial~correlation} \quad $'+windowtxt)
      ax.set_xlabel(xlabel)

      itheta=0

      y = np.sum(self.phisq[:,itheta,1:,:],axis=1)
      for j in range(self.n_radial):
         ave[j] = average(y[j,:],self.t,w)

      ave = np.roll(ave,-self.n_radial/2)
      ave[0] = 0.0
      corr = np.fft.fft(ave,self.n_radial)
      corr = np.fft.fftshift(corr)
      corr /= np.max(np.abs(corr))
      corr = corr.real
      delta_r = np.fft.fftfreq(self.n_radial)
      delta_r = np.fft.fftshift(delta_r)
      Lx = 2*np.pi/dk
      delta_r *= Lx

      # calculate envelope
      corr_hilbert = signal.hilbert(corr)
      corr_env = np.abs(corr_hilbert)
      ax.set_ylabel(r'$C_{\delta \phi}(\Delta r)$',color='k')
      ax.plot(delta_r, 0*delta_r, color='k', ls='--')
      ax.plot(delta_r, corr, color=color[0])

      l_corr, pcov = curve_fit(absexp, delta_r, corr_env, p0=10.0)
      ax.plot(delta_r,absexp(delta_r,l_corr),color=color[1],ls='-.')

      ax.set_xlim([np.min(delta_r),np.max(delta_r)])
      ax.set_ylim(-1,1)

      print 'INFO: (data_plot.py) l_corr = ',l_corr


   def plot_geo(self,fig=None):

      self.getgeo()

      if fig is None:
         fig = plt.figure(figsize=(15,9))
         fig.subplots_adjust(left=0.07,right=0.95,top=0.94,bottom=0.06,wspace=0.25,hspace=0.32)

      theta = self.geo[:,0]/np.pi

      for p in range(8):
         p1 = p+1
         y = self.geo[:,p1]

         ax = fig.add_subplot(2,4,p1)
         ax.grid(which="majorminor",ls=":")
         ax.grid(which="major",ls=":")
         ax.set_xlabel(r'$\theta/\pi$')
         ax.set_title(r'$'+self.geotag[p1]+'$')
         ax.plot(theta,y)
         ax.set_xlim([-1,1])
