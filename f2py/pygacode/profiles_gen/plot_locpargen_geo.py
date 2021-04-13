import os
import sys
import wx
import matplotlib
import numpy as np
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import matplotlib.pyplot as plt

matplotlib.rc('text',usetex=True)
matplotlib.rc('font',size=18)
ntheta = int(sys.argv[1])
print(ntheta)

wdir = os.path.realpath('./')

data = np.fromfile('bin.locpargen.theta',dtype='float32')
data = np.reshape(data,(6,ntheta),'F')
x = data[0,:]/np.pi

def plot_select(ax,tag):

    ax.set_xlabel(r'$\theta$')
    ax.grid(which="both",ls=":")
    ax.grid(which="major",ls=":")
    
    if tag == 'Jr':
        rr=data[1,:]
        rt=data[2,:]
        zr=data[3,:]
        zt=data[4,:]

        ax.plot(x,rr*zt-zr*rt,color='k',label=r'$R_r Z_\theta - Z_r R_\theta$')
        ax.plot(x,rr*zt,color='r',alpha=0.4,label=r'$R_r Z_\theta$')
        ax.plot(x,-zr*rt,color='b',alpha=0.4,label=r'$-Z_r R_\theta$')
        ax.plot(x,x*0.0,linestyle='--',color='k')
        ax.legend()
    elif tag == 'gsin':
        rr=data[:,1]
        rt=data[:,2]
        zr=data[:,3]
        zt=data[:,4]

        ax.plot(x,rr*zt-zr*rt,color='k',label=r'$R_r Z_\theta - Z_r R_\theta$')
        ax.plot(x,rr*zt,color='r',alpha=0.4,label=r'$R_r Z_\theta$')
        ax.plot(x,-zr*rt,color='b',alpha=0.4,label=r'$-Z_r R_\theta$')
        ax.plot(x,x*0.0,linestyle='--')
        ax.legend()

        
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

        mytabs = ['Jr','gsin']

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
