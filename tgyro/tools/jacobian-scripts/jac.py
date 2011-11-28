#!/usr/bin/python
#!/usr/bin/python
#try to color some plots
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.text as text
import matplotlib.mlab as mlab
import matplotlib.cbook as cbook

def plot_jac():
  
    data = np.genfromtxt('output/jac.out')
    vecdata = np.genfromtxt('output/vec.out')
    
    dqidzi=[]
    dqidze=[]
    dqedzi=[]
    dqedze=[]

    z_i=[]
    q_i=[]
    z_e=[]
    q_e=[]
    
    for i in range(0,len(data),4):
        dqidzi.append(data[i])
        dqidze.append(data[i+1])
        dqedzi.append(data[i+2])
        dqedze.append(data[i+3])
                      

    for i in range(0,4):
        z_i.append(vecdata[i])
        q_i.append(vecdata[i+4])
        z_e.append(vecdata[i+8])
        q_e.append(vecdata[i+12])

        
    lndqidzi=[]
    lndqidze=[]
    lndqedzi=[]
    lndqedze=[]

    for i in range(0,4):
        lndqidzi.append(1/q_i[i]*dqidzi[i])
        lndqidze.append(1/q_i[i]*dqidze[i])
        lndqedzi.append(1/q_e[i]*dqedzi[i])
        lndqedze.append(1/q_e[i]*dqedze[i])


    lndqidlnzi=[]
    lndqidlnze=[]
    lndqedlnzi=[]
    lndqedlnze=[]

    for i in range(0,4):
        lndqidlnzi.append(1/q_i[i]*dqidzi[i]*z_i[i])
        lndqidlnze.append(1/q_i[i]*dqidze[i]*z_e[i])
        lndqedlnzi.append(1/q_e[i]*dqedzi[i]*z_i[i])
        lndqedlnze.append(1/q_e[i]*dqedze[i]*z_e[i])


    r0=0.1875
    
    r=[]

    for i in range(0,4):
        r.append(r0*(i+1))
        
#    print 'radius', r
 


    prop = matplotlib.font_manager.FontProperties(size=14)
    fig = plt.figure(figsize=(6,4))
    plot = fig.add_subplot(111)
    plot2 = fig.add_subplot(111)



    plot.set_xlabel(r'r/a')
    plot.set_ylabel(r'$dQ_\sigma/dz_\alpha$')
    plt.title('Jacobian')
    plot.grid(True)
    
    plot.plot(r,dqidzi,'-o',color='darkgreen', lw=2.5, label="dQidzi")
    plot.plot(r,dqidze,'-o',color='blue', lw=2.5, label="dQidze")
    plot.plot(r,dqedzi,'-o',color='red', lw=2.5, label="dQedzi")
    plot.plot(r,dqedze,'-o',color='magenta', lw=2.5, label="dQedze")

    handles, labels=plot.get_legend_handles_labels()
    plot.legend(labels, 'best')

    plt.savefig('plots/dQdz.pdf')




if __name__=="__main__":
    print("Plotting jacobian")
    plot_jac()

