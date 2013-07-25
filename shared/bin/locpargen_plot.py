import wx
import matplotlib
import numpy as np
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from gacodeplotdefs import *

data = np.loadtxt('out.locpargen.geo')
data2 = np.loadtxt('out.locpargen.geo.2')

class TabPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent=parent)
 
        self.figure = plt.Figure()
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()

    def draw(self,i,ylabel):

        self.ax.grid(which="majorminor",ls=":")
        self.ax.grid(which="major",ls=":")
        self.ax.set_xlabel(r'$\theta/\pi$',fontsize=GFONTSIZE)
        self.ax.set_ylabel(r'$'+ylabel+'$',color='k',fontsize=GFONTSIZE)
        #self.ax.plot(data[:,0]/np.pi,data[:,i])
        self.ax.plot(data2[:,0]/np.pi,data2[:,i],color='m')
        self.ax.set_xlim([-1,1])
        if ylabel == '\\Theta':
            self.ax.plot(data[:,0]/np.pi,data[:,0],'--k')
           
            
class DemoFrame(wx.Frame):
     def __init__(self):
        """Constructor"""        
        wx.Frame.__init__(self, None, wx.ID_ANY, 
                          "Geometry functions",
                          size=(1000,800))
        panel = wx.Panel(self)
 
        notebook = wx.Notebook(panel)

        tab1 = TabPanel(notebook)
        tab1.draw(1,'B')
        notebook.AddPage(tab1,"b")
 
        tab2 = TabPanel(notebook)
        tab2.draw(2,'\\partial B/\\partial \\theta')
        notebook.AddPage(tab2,"dbdt")
 
        tab3 = TabPanel(notebook)
        tab3.draw(3,'\\partial^2 B / \\partial \\theta^2')
        notebook.AddPage(tab3,"dbdt2")

        tab4 = TabPanel(notebook)
        tab4.draw(4,'B_p')
        notebook.AddPage(tab4,"bp")

        tab5 = TabPanel(notebook)
        tab5.draw(5,'R')
        notebook.AddPage(tab5,"bigr")

        tab6 = TabPanel(notebook)
        tab6.draw(6,'| \\nabla R |')
        notebook.AddPage(tab6,"grad_r")

        tab7 = TabPanel(notebook)
        tab7.draw(7,'J_r')
        notebook.AddPage(tab7,"jac_r")

        tab8 = TabPanel(notebook)
        tab8.draw(8,'\\Theta')
        notebook.AddPage(tab8,"captheta")
 
        tab9 = TabPanel(notebook)
        tab9.draw(9,'\\partial \ell / \\partial \\Theta')
        notebook.AddPage(tab9,"l_t")
 
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(notebook, 1, wx.ALL|wx.EXPAND, 5)
        panel.SetSizer(sizer)
        self.Layout()
 
        self.Show()
 
if __name__ == "__main__":
    app = wx.App(False)
    frame = DemoFrame()
    app.MainLoop()
