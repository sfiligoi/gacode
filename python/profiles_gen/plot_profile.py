#-------------------------------------------------------------
# plot_profile.py
#
# PURPOSE:
#  Notebook plotter to see input.gacode profiles
#-------------------------------------------------------------

import os
import wx
import sys
import matplotlib
import numpy as np
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import matplotlib.pyplot as plt
from gacodefuncs import *
from gacode import expro

matplotlib.rc('text',usetex=True)
matplotlib.rc('font',size=18)

simdir = './'
wdir = os.path.realpath(simdir)
rvar = sys.argv[1]
rmin = sys.argv[2]
rmax = sys.argv[3]
ext = sys.argv[4]
loc = int(sys.argv[5])

m1 = 0 ; m2 = 0

def plotit(ax,x,y,ystr):
   global m1,m2
   ax.plot(x[m1:m2],y[m1:m2],label=r'$'+ystr+'$')
   return

def plot_select(ax,tag):
   global m1,m2
   
   expro.expro_read('input.gacode')

   x = expro.expro_rmin ; x = x/max(x)
   n = expro.expro_n_ion
   
   # normalization
   csa = expro.expro_cs/expro.expro_rmin[-1]

   sname = gapystr(expro.expro_name)
   stype = gapystr(expro.expro_type)

   # Set x-range
   m1 = 0 ; m2 = len(x)
   if rmax != 'auto':
      m2 = np.argmin(np.abs(x-np.float(rmax)))
   if rmin != 'auto':
      m1 = np.argmin(np.abs(x-np.float(rmin)))
      
   ax.set_xlim([x[m1],x[m2-1]])

   ax.set_xlabel(r'$r/a$')

   m=m2
   if tag == 'bunit':
      # bunit
      y = expro.expro_bunit ; ystr = 'B_\mathrm{unit}' ; plotit(ax,x,y,ystr)

   if tag == 'gammae':
      # gamma_e
      y = expro.expro_gamma_e ; ystr = '(a/c_s) \gamma_E' ; plotit(ax,x,y/csa,ystr)

   if tag == 'gammap':
      # gamma_p
      y = expro.expro_gamma_p ; ystr = '(a/c_s) \gamma_p' ; plotit(ax,x,y/csa,ystr) 

   if tag == 'mach':
      # mach
      y = expro.expro_mach ; ystr = r'M' ; plotit(ax,x,y,ystr)

   if tag == 'r':
      # rho
      y = expro.expro_rho ; ystr = '\\rho' ; plotit(ax,x,y,ystr)
      # polflux
      y = expro.expro_polflux ; ystr = '\psi' ; plotit(ax,x,y/y[-1],ystr)

   if tag == 'n':
      # ne
      y = expro.expro_ne ; ystr = 'n_e' ; plotit(ax,x,y,ystr)
      # ni
      for p in range(n):
         y = expro.expro_ni[p,:] ; ystr = 'n_i ['+sname[p]+']' ; plotit(ax,x,y,ystr)

   if tag == 'T':
      # Te
      y = expro.expro_te ; ystr = 'T_e' ; plotit(ax,x,y,ystr)
      # Ti
      for p in range(n):
         y = expro.expro_ti[p,:] ; ystr = 'T_i ['+sname[p]+']' ; plotit(ax,x,y,ystr)
         
   if tag == 'kappa':
       y = expro.expro_kappa ; ystr = '\kappa' ; plotit(ax,x,y,ystr)
       y = expro.expro_skappa ; ystr = 's_\kappa' ; plotit(ax,x,y,ystr)

   if tag == 'sin':
       y = np.arcsin(expro.expro_delta) ; ystr = 's_1 = \sin^{-1}\delta' ; plotit(ax,x,y,ystr)
       y = -expro.expro_zeta         ; ystr = 's_2 = -\zeta' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_sin3    ; ystr = 's_3' ; plotit(ax,x,y,ystr)
     
   if tag == 'cos':
       y = expro.expro_shape_cos0 ; ystr = 'c_0' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_cos1 ; ystr = 'c_1' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_cos2 ; ystr = 'c_2' ; plotit(ax,x,y,ystr)
       y = expro.expro_shape_cos3 ; ystr = 'c_3' ; plotit(ax,x,y,ystr)
     
   if tag == 'q':
       y = expro.expro_q ; ystr = 'q'
       if y[0] < 0.0:
          y = -y
          ystr = '-'+ystr
       plotit(ax,x,y,ystr)
       y = expro.expro_s ; ystr = 's' ; plotit(ax,x,y,ystr)

   ax.legend(loc=loc)
       

#-------------------------------------------------------------------------------------

class TabPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent=parent)
 
        self.figure = plt.Figure()
        self.figure.subplots_adjust(left=0.07,right=0.95)
        self.ax = self.figure.add_subplot(111)
        self.ax.grid(which="both",ls=":")
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
                          size=(1100,800))

        panel = wx.Panel(self)
 
        notebook = wx.Notebook(panel)

        tab = TabPanel(notebook)
        tab.draw('gammae')
        notebook.AddPage(tab,'gammae')

        tab = TabPanel(notebook)
        tab.draw('gammap')
        notebook.AddPage(tab,'gammap')

        tab = TabPanel(notebook)
        tab.draw('mach')
        notebook.AddPage(tab,'mach')

        tab = TabPanel(notebook)
        tab.draw('r')
        notebook.AddPage(tab,'r')

        tab = TabPanel(notebook)
        tab.draw('kappa')
        notebook.AddPage(tab,'kappa')


        tab = TabPanel(notebook)
        tab.draw('bunit')
        notebook.AddPage(tab,'bunit')

        tab = TabPanel(notebook)
        tab.draw('q')
        notebook.AddPage(tab,'q')

        tab = TabPanel(notebook)
        tab.draw('n')
        notebook.AddPage(tab,'n')

        tab = TabPanel(notebook)
        tab.draw('T')
        notebook.AddPage(tab,'T')

        tab = TabPanel(notebook)
        tab.draw('sin')
        notebook.AddPage(tab,'sin')
        tab = TabPanel(notebook)
        tab.draw('cos')
        notebook.AddPage(tab,'cos')

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
