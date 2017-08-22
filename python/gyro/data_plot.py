import data
import sys
import numpy as np
from gacodeplotdefs import *
from gyro.data import GYROData

#---------------------------------------------------------------------------#
class gyrodata_plot(data.GYROData):
    def plot_freq(self, sim='.', window=0.5, fig=None):
        """Plot frequency vs time"""
        """window = avg. window fraction"""
        
        GFONTSIZE=16
        if (fig is None):
            fig = plt.figure()
        
        t    = sim.t['(c_s/a)t']

        # Read freq data
        sim.read_freq()

        # Determine tmin
        for i in range(len(t)):
            if t[i] < (1.0-window)*t[len(t)-1]:
                imin = i

        color = ['k','m','b','c']
        tor_n = self.n0 + self.dn*arange(0,self.n_n)
        #======================================
        ax = fig.add_subplot(121)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(r'$(c_s/a) t$',fontsize=GFONTSIZE)
        ax.set_ylabel(r'$(a/c_s)\gamma$',color='k',fontsize=GFONTSIZE)
        #=====================================

        #Gamma 
        for i in range(sim.profile['n_n']):
            ax.plot(t[imin:],sim.freq['(a/c_s)gamma'][i,imin:],color=color[i],
                    label='gamma_n%d'%tor_n[i])

        #======================================
        ax = fig.add_subplot(122)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(r'$(c_s/a) t$',fontsize=GFONTSIZE)
        ax.set_ylabel(r'$(a/c_s)\omega$',color='k',fontsize=GFONTSIZE)
        #=====================================

        #Omega
        for i in range(sim.profile['n_n']):
            ax.plot(t[imin:],sim.freq['(a/c_s)w'][i,imin:],color=color[i],
                    label='omega_n%d'%tor_n[i])

