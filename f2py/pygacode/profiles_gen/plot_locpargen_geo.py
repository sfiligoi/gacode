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

wdir = os.path.realpath('./')

data = np.fromfile('bin.locpargen.theta',dtype='float32')
nvec = len(data)//ntheta
data = np.reshape(data,(ntheta,nvec),'F')
x = data[:,0]/np.pi

def plot_select(ax,tag):

    ax.set_xlabel(r'$\theta/\pi$')
    ax.grid(which="both",ls=":")
    ax.grid(which="major",ls=":")

    ax.plot(x,x*0.0,linestyle='--',color='k')
    
    if tag == 'gradr':
        y=data[:,1]
        ax.plot(x,y,color='k',label=r'$\left| \nabla r \right|$')

    elif tag == 'R':
        y=data[:,2]
        ax.plot(x,y,color='k',label=r'$R$')

    elif tag == 'Z':
        y=data[:,3]
        ax.plot(x,y,color='k',label=r'$Z$')

    elif tag == 'gsin':
        y=data[:,4]
        ax.plot(x,y,color='k',label=r'$\mathrm{gsin}$')

    elif tag == 'gcos1':
        y=data[:,5]
        ax.plot(x,y,color='k',label=r'$\mathrm{gcos}_1$')

    elif tag == 'gcos2':
        y=data[:,6]
        ax.plot(x,y,color='k',label=r'$\mathrm{gcos}_2$')
        
    elif tag == 'Bt':
        y=data[:,7]
        ax.plot(x,y,color='k',label=r'$B_t$')
 
    elif tag == 'Bp':
        y=data[:,8]
        ax.plot(x,y,color='k',label=r'$B_p$')

    elif tag == 'g_theta':
        y=data[:,9]
        ax.plot(x,y,color='k',label=r'$G_\theta$')

    elif tag == 'gq':
        y=data[:,10]
        ax.plot(x,y,color='k',label=r'$G_q$')

    elif tag == 'captheta':
        y=data[:,11]
        ax.plot(x,y,color='k',label=r'$\Theta$')

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

        mytabs = ['gradr','R','Z','gsin','gcos1','gcos2','Bt','Bp','g_theta','gq','captheta']

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
