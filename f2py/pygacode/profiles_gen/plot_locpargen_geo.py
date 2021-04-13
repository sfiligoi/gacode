import os
import wx
import matplotlib
import numpy as np
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import matplotlib.pyplot as plt
#from ..gacodefuncs import *

matplotlib.rc('text',usetex=True)
matplotlib.rc('font',size=18)

data = np.genfromtxt('out.locpargen.theta')
simdir = './'
wdir = os.path.realpath(simdir)

def plot_select(ax,tag):
    
    if tag == 'J':
        x =data[:,0]/np.pi
        rr=data[:,1]
        rt=data[:,2]
        zr=data[:,3]
        zt=data[:,4]

        ax.set_xlabel(r'$\theta$')
        ax.set_ylabel(r'$F$')
        ax.grid(which="both",ls=":")
        ax.grid(which="major",ls=":")

        ax.plot(x,rr*zt-zr*rt,color='k')
        ax.plot(x,rr*zt,color='r')
        ax.plot(x,-zr*rt,color='b')
        ax.plot(x,x*0.0,linestyle='--')

    ax.set_xlim([-1,1])
    plt.tight_layout()
    
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
                          'locpargen geo plot -- '+wdir,
                          size=(1100,800))

        panel = wx.Panel(self)
 
        notebook = wx.Notebook(panel)

        mytabs = ['J']

        for x in mytabs:
           tab = TabPanel(notebook)
           tab.draw(x)
           notebook.AddPage(tab,x)

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
