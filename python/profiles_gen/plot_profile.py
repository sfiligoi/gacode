#-------------------------------------------------------------
# plot_profile.py
#
# PURPOSE:
#  Notebook plotter to see tgyro results.
#-------------------------------------------------------------

import os
import wx
import matplotlib
import sys
import numpy as np
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib import rc
from gacodefuncs import *
from tgyro.data import tgyrodata
from matplotlib.colors import LogNorm
from gapy import expro

rc('text',usetex=True)
rc('font',size=18)

simdir = './'
wdir = os.path.realpath(simdir)
  
def gapystr(s):

   n = len(s[0])
   u = s.transpose().reshape(-1,n).view('S'+str(n))
   
   return [u[0].tostring().strip(),
           u[1].tostring().strip(),
           u[2].tostring().strip(),
           u[3].tostring().strip()]

def plot_select(ax,tag):

   # Helper routine to plot data (tag) from input.profiles

   expro.expro_read('input.gacode')

   x = expro.expro_rmin ; x = x/max(x)
   n = expro.expro_n_ion
   
   if tag == 'pro':
      y = expro.expro_ne ; ystr = expro.expro_ne_str[2].strip()
      ax.plot(x,y,label=r'$'+ystr+'$')
      for p in range(n):
         y = expro.expro_ni ; ystr = expro.expro_ni_str[2].strip()
         ax.plot(x,y[p,:],label=r'$'+ystr+str(p+1)+'}$')
      ax.legend()
         
   if tag == 'geo':
       y = expro.expro_rmaj ; ystr = expro.expro_rmaj_str[2].strip()
       ax.plot(x,y,label=r'$'+ystr+'$')

       y = expro.expro_zmag ; ystr = expro.expro_zmag_str[2].strip()
       ax.plot(x,y,label=r'$'+ystr+'$')

       y = expro.expro_kappa ; ystr = expro.expro_kappa_str[2].strip()
       ax.plot(x,y,label=r'$'+ystr+'$')
     
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
        tab.draw('pro')
        notebook.AddPage(tab,'pro')

        tab = TabPanel(notebook)
        tab.draw('geo')
        notebook.AddPage(tab,'geo')

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
