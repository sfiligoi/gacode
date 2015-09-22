#-------------------------------------------------------------
# plot_notebook.py
#
# PURPOSE:
#  Notebook of major cgyro plots.
#-------------------------------------------------------------

import wx
import matplotlib
import string
import sys
import re
import numpy as np
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from gacodeplotdefs import *
from cgyro.data import cgyrodata

w     = float(sys.argv[1])
ext   = sys.argv[2]
is0   = int(sys.argv[3])

sim = cgyrodata('./')
sim.getbig()

ns = sim.n_species
t = sim.t
 
# Determine tmin
for i in range(len(t)):
    if t[i] < (1.0-w)*t[len(t)-1]:
        imin = i

tstr = r'$'+str(sim.t[imin])+' < (c_s/a) t < '+str(sim.t[-1])+'$' 
color = ['m','k','b','c']

#-------------------------------------------------------------------------------------

def plot_gen(ax,tag):
    
    if tag == 'fluxqt':

        y = np.sum(sim.flux_e,axis=(0,2))
        ax.set_title(tstr)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(r'$(c_s/a)\, t$')
        for ispec in range(ns):
            ave = average(y[ispec,:],t,w)
            y_ave = ave*np.ones(len(t))
            label = r'$Q_'+str(ispec)+': '+str(round(ave,3))+'$'
            # Average
            ax.plot(t[imin:],y_ave[imin:],'--',color=color[ispec])
            # Time trace
            ax.plot(sim.t,y[ispec,:],label=label,color=color[ispec])

        ax.legend(loc=2)

    elif tag == 'fluxgt':
        
        y = np.sum(sim.flux_n,axis=(0,2))
        ax.set_title(tstr)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(r'$(c_s/a)\, t$')
        ax.set_title(tstr)
        for ispec in range(ns):
            ave = average(y[ispec,:],t,w)
            y_ave = ave*np.ones(len(t))
            label = r'$Q_'+str(ispec)+': '+str(round(ave,3))+'$'
            # Average
            ax.plot(t[imin:],y_ave[imin:],'--',color=color[ispec])
            # Time trace
            ax.plot(sim.t,y[ispec,:],label=label,color=color[ispec])

        ax.legend(loc=2)
 
    elif tag == 'phinp':

        nx=sim.n_radial
        ny=sim.n_n
        f = np.zeros([nx,ny])
        n = sim.n_time

        for i in np.arange(imin,n):
            f = f+sim.phisq[:,:,i]

        f = 1e-12+f/(n-imin)
        f = np.log10(f)
            
        ax.set_xlabel(r'$k_y \rho_s$',fontsize=GFONTSIZE)
        ax.set_ylabel(r'$k_x \rho_s$',fontsize=GFONTSIZE)
        ax.set_title(tstr)
        fmax = f.max()
        fmin = f.max()-5

        d = (fmax-fmin)/200.0
        levels = np.arange(fmin-d,fmax+d,d)
            
        ax.contourf(sim.ky,sim.kx,f,levels,origin='lower')

    elif tag == 'nnp':

        nx=sim.n_radial
        ny=sim.n_n
        f = np.zeros([nx,ny])
        n = sim.n_time

        for i in np.arange(imin,n):
            f = f+sim.nsq[:,is0,:,i]

        f = 1e-12+f/(n-imin)
        f = np.log10(f)
            
        ax.set_xlabel(r'$k_y \rho_s$',fontsize=GFONTSIZE)
        ax.set_ylabel(r'$k_x \rho_s$',fontsize=GFONTSIZE)
        ax.set_title(tstr)
        fmax = f.max()
        fmin = f.max()-5

        d = (fmax-fmin)/200.0
        levels = np.arange(fmin-d,fmax+d,d)
            
        ax.contourf(sim.ky,sim.kx,f,levels,origin='lower')

    elif tag == 'fluxqnpt':

        nx=sim.n_radial
        ny=sim.n_n
        f = np.zeros([nx,ny])
        n = sim.n_time

        for i in np.arange(imin,n):
            f = f+sim.flux_e[:,is0,:,i]

        f = f/(n-imin)
            
        ax.set_xlabel(r'$k_y \rho_s$',fontsize=GFONTSIZE)
        ax.set_ylabel(r'$k_x \rho_s$',fontsize=GFONTSIZE)
        ax.set_title(tstr)
        fmax = f.max()
        fmin = f.min()

        d = (fmax-fmin)/200.0
        levels = np.arange(fmin-d,fmax+d,d)
            
        ax.contourf(sim.ky,sim.kx,f,levels,origin='lower')

        
#-------------------------------------------------------------------------------------

class TabPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent=parent)
 
        self.figure = plt.Figure()
        self.figure.subplots_adjust(left=0.07,right=0.95)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()

    def draw(self,tag):
        plot_gen(self.ax,tag)
        
#-------------------------------------------------------------------------------------

class DemoFrame(wx.Frame):
  def __init__(self):
        """Constructor"""        
        wx.Frame.__init__(self, None, wx.ID_ANY, 
                          "CGYRO Simulation Summary",
                          size=(1100,600))
        panel = wx.Panel(self)
 
        notebook = wx.Notebook(panel)

        tab = TabPanel(notebook)
        tab.draw('fluxqt')
        notebook.AddPage(tab,'Q(t)')

        tab = TabPanel(notebook)
        tab.draw('fluxqnpt')
        notebook.AddPage(tab,'Q_np(t)')

        tab = TabPanel(notebook)
        tab.draw('fluxgt')
        notebook.AddPage(tab,'Gamma(t)')

        tab = TabPanel(notebook)
        tab.draw('phinp')
        notebook.AddPage(tab,'|Phi|_np')

        tab = TabPanel(notebook)
        tab.draw('nnp')
        notebook.AddPage(tab,'|n|_np')

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(notebook, 1, wx.ALL|wx.EXPAND, 5)
        panel.SetSizer(sizer)
        self.Layout()
 
        self.Show()
 
#-------------------------------------------------------------------------------------

if __name__ == "__main__":

  if ext == 'screen':

    # On-screen wxpython notebook
      
    app = wx.App(False)
    frame = DemoFrame()
    app.MainLoop()
    
  else:

    # Generate plots

    list=['fluxqt','fluxgt','phinp','nnp','fluxqnpt']

    for x in list:
        figure = plt.figure(figsize=(8,5))
        figure.subplots_adjust(left=0.11,right=0.95,bottom=0.16)
        ax = figure.add_subplot(111)
        plot_gen(ax,x)
        pfile = 'out.'+x+'.'+ext
        plt.savefig(pfile)
        print 'INFO: (plot_notebook.py) Wrote '+pfile
