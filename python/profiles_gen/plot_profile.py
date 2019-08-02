#-------------------------------------------------------------
# plot_profile.py
#
# PURPOSE:
#  Notebook plotter to see input.gacode profiles
#-------------------------------------------------------------

import os
import wx
import matplotlib
import sys
import numpy as np
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import matplotlib.pyplot as plt
from gacodefuncs import *
from pygacode import expro

matplotlib.rc('text',usetex=True)
matplotlib.rc('font',size=18)

simdir = './'
wdir = os.path.realpath(simdir)
rvar = sys.argv[1]
rmin = sys.argv[2]
rmax = sys.argv[3]
ext = sys.argv[4]
loc = int(sys.argv[5])

def gapystr(s):

   if sys.version_info[0] > 2:
      # Python 3 
      return s[2].decode('UTF-8')
   else:
      # Python 2
      n = len(s[0])
      u = s.transpose().reshape(-1,n).view('S'+str(n))   
      return u[2].tostring().strip()

def gapystrv(s):

   if sys.version_info[0] > 2:
      # Python 3
      u=[]
      for i in range(expro.expro_n_ion):
         u.append(s[i].decode('UTF-8').strip())
      return u
   else:
      # Python 2
      n = len(s[0])
      u = s.transpose().reshape(-1,n).view('S'+str(n))   
      return u[:].tostring().strip()
      
def plot_select(ax,tag):

   # Helper routine to plot data (tag) from input.gacode

   expro.expro_read('input.gacode')

   x = expro.expro_rmin ; x = x/max(x)
   n = expro.expro_n_ion
   
   # normalization
   csa = expro.expro_cs/expro.expro_rmin[-1]

   # Set x-range
   m = len(x)
   if rmax != 'auto':
      ax.set_xlim([0,float(rmax)])
      for m in range(len(x)):
         if x[m] > float(rmax):
            break

   ax.set_xlabel(r'$r/a$')
      
   if tag == 'bunit':
      # bunit
      y = expro.expro_bunit ; ystr = r'$B_\mathrm{unit}$'
      ax.plot(x[:m],y[:m],label=r'$'+ystr+'$')

   if tag == 'gammae':
      # gamma_e
      y = expro.expro_gamma_e ; ystr = gapystr(expro.expro_gamma_e_str)
      ax.plot(x[:m],y[:m]/csa[:m],label=r'$'+ystr+'$')

   if tag == 'gammap':
      # gamma_p
      y = expro.expro_gamma_p ; ystr = gapystr(expro.expro_gamma_p_str)
      ax.plot(x[:m],y[:m]/csa[:m],label=r'$'+ystr+'$')

   if tag == 'r':
      # rho
      y = expro.expro_rho ; ystr = gapystr(expro.expro_rho_str)
      ax.plot(x[:m],y[:m],label=r'$'+ystr+'$')
      # polflux
      y = expro.expro_polflux ; ystr = gapystr(expro.expro_polflux_str)
      ax.plot(x[:m],y[:m]/y[-1],label=r'$'+ystr+'$')

   if tag == 'n':
      sname = gapystrv(expro.expro_name)
      stype = gapystrv(expro.expro_type)
      # ne
      y = expro.expro_ne ; ystr = gapystr(expro.expro_ne_str)
      ax.plot(x[:m],y[:m],label=r'$'+ystr+'$')
      # ni
      y = expro.expro_ni ; ystr = gapystr(expro.expro_ni_str)
      for p in range(n):
         ax.plot(x[:m],y[p,:m],label=r'$'+ystr+sname[p]+'}$')

   if tag == 'T':
      sname = gapystrv(expro.expro_name)
      stype = gapystrv(expro.expro_type)
      # Te
      y = expro.expro_te ; ystr = gapystr(expro.expro_te_str)
      ax.plot(x[:m],y[:m],label=r'$'+ystr+'$')
      # Ti
      y = expro.expro_ti ; ystr = gapystr(expro.expro_ti_str)
      for p in range(n):
         if stype[p] == '[therm]':
            ax.plot(x[:m],y[p,:m],label=r'$'+ystr+sname[p]+'}$')
         
   if tag == 'kappa':
       y = expro.expro_kappa ; ystr = gapystr(expro.expro_kappa_str)
       ax.plot(x[:m],y[:m],label=r'$'+ystr+'$')
       y = expro.expro_skappa ; ystr = 's_\kappa'
       ax.plot(x[:m],y[:m],label=r'$'+ystr+'$')

   if tag == 'delta':
       y = expro.expro_delta ; ystr = gapystr(expro.expro_delta_str)
       ax.plot(x[:m],y[:m],label=r'$'+ystr+'$')
       y = expro.expro_sdelta ; ystr = 's_\delta'
       ax.plot(x[:m],y[:m],label=r'$'+ystr+'$')

   ax.legend()
       

      
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
        plot_select(self.ax,tag)
  
#-------------------------------------------------------------------------------------

class DemoFrame(wx.Frame):
   
  def __init__(self):

        wx.Frame.__init__(self, None, wx.ID_ANY, 
                          'input.gacode plotting notebook -- '+wdir,
                          size=(1100,600))

        panel = wx.Panel(self)
 
        notebook = wx.Notebook(panel)

        tab = TabPanel(notebook)
        tab.draw('bunit')
        notebook.AddPage(tab,'bunit')

        tab = TabPanel(notebook)
        tab.draw('gammae')
        notebook.AddPage(tab,'gammae')

        tab = TabPanel(notebook)
        tab.draw('gammap')
        notebook.AddPage(tab,'gammap')

        tab = TabPanel(notebook)
        tab.draw('r')
        notebook.AddPage(tab,'r')

        tab = TabPanel(notebook)
        tab.draw('kappa')
        notebook.AddPage(tab,'kappa')

        tab = TabPanel(notebook)
        tab.draw('delta')
        notebook.AddPage(tab,'delta')

        tab = TabPanel(notebook)
        tab.draw('n')
        notebook.AddPage(tab,'n')

        tab = TabPanel(notebook)
        tab.draw('T')
        notebook.AddPage(tab,'T')

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(notebook, 1, wx.ALL|wx.EXPAND, 5)
        panel.SetSizer(sizer)
        self.Layout() 
        self.Show()
 

if __name__ == "__main__":

   # On-screen wxpython notebook
      
   app = wx.App(False)
   frame = DemoFrame()
   app.MainLoop()
