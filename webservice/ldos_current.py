def get_band(lead0):
    import numpy as np
    import kwant
    import matplotlib.pyplot as plt
    bands = kwant.physics.Bands(lead0.finalized())
    t=2.75
    momenta = np.linspace(-np.pi, np.pi, 101)
    energies = [bands(k) for k in momenta]
    for j in range(len(energies)):
        energies[j]=energies[j]/t  
    for j in range(len(momenta)):
        momenta[j]=momenta[j]/np.pi 
    plt.figure(figsize=(10,10))
    plt.plot(momenta, energies)
    plt.xticks(np.arange(-1,1.1,1), ('X', 'Γ', 'X'))
    plt.xlabel("Wave vector")
    plt.ylabel("Energy (t)")
    plt.grid(b=None, which='major', axis='both')
    plt.savefig('webservice/user_static/img/band.png',bbox_inches='tight')

def get_ldos(fsys,xaxis,yaxis,en=0.5,t=2.75):
    import kwant
    import matplotlib.pyplot as plt
    import numpy as np
    ldos=kwant.ldos(fsys,t*(en+0.0001))
    max_ldos=max(ldos)
    ind = np.argmax(ldos)
    ldos=ldos/max(ldos)
    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_xticks([])
    ax.set_yticks([])    
    kwant.plot(fsys,site_size=ldos/1.75,site_color='blue', lead_color='black',ax=ax)
    circle1 = plt.Circle((xaxis[0]+0.5,yaxis[1]-0.5), 0.325, color='blue')
    ax.add_artist(circle1)
    ax.text(xaxis[0]+1,yaxis[1]-0.6,s=str("{:.2E}".format(max_ldos)))
    ax.set_aspect('equal', 'box')
    ax.grid(b=None, which='major', axis='both')
    fig.savefig('webservice/user_static/img/ldos.png',bbox_inches='tight')

def lead_wfn(perfsys,en=0.5,mode=0,t=2.75):
    import kwant
    import numpy as np
    from cmath import phase
    from math import pi
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import cmocean
    import unicodedata
    from matplotlib import gridspec
    leads=perfsys
    ######Calculating the wave function in scattering center
    wf = kwant.solvers.default.wave_function(leads,energy=t*(en+0.00001))
    psi_l= wf(0)[mode] #selecting the wfn incoming from left (0) and selecting the mode
    amplitude=[]
    pha=[]
    for j in range (len(psi_l)): #for each atom get amplitude and phase
        amplitude.append((np.abs(psi_l[j])))
        pha.append(phase(psi_l[j])/pi)
    first=0 
    for count,amp in enumerate(amplitude): #add normalisation for the correct color map range (No idea for a better solution)
        if amp==0.0: #select the "fake" atoms and assign a phase of -1 or 1
            if first==1:
                pha[count]=-1
            if first==0:
                first=1
                pha[count]=1
    fig = plt.figure(figsize=(10,3)) # create the canvas for plotting
    gs = gridspec.GridSpec(1, 2, width_ratios=[9, 1]) 
    ax1 = plt.subplot(gs[0]) 
    ax2 = plt.subplot(gs[1],projection='polar')    
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_aspect('equal','box')
    kwant.plotter.plot(leads,site_size=amplitude/max(amplitude)/2.5,site_color=pha,lead_color='black',cmap=cmocean.cm.phase,fig_size=(10,5),colorbar=False,ax=ax1)
    azimuths = np.arange(0, 361, 1)
    zeniths = np.arange(70, 100, 1)
    values = azimuths * np.ones((30, 361))    
    ax2.pcolormesh(azimuths*np.pi/180.0, zeniths, values, cmap=cmocean.cm.phase)
    ax2.set_xticklabels(['0', '', unicodedata.lookup("GREEK SMALL LETTER PI")+"/2", '', unicodedata.lookup("GREEK SMALL LETTER PI"), '', "-"+unicodedata.lookup("GREEK SMALL LETTER PI")+"/2", ''])
    ax2.set_yticks([])    
    fig.savefig('webservice/user_static/img/wfn.png',bbox_inches='tight' )                                                                                                                                                                                                
def get_current(fsys,xaxis,yaxis,en=0.5,t=2.75):
    import kwant
    import numpy as np
    import matplotlib as mpl
    from matplotlib import gridspec
    import matplotlib.pyplot as plt
    import operator
    from matplotlib import cm
    psi = kwant.wave_function(fsys,energy=t*(en+0.001))(0)
    hop_lead_nr=kwant.plotter.sys_leads_hoppings(fsys,num_lead_cells=0)
    hoppings=kwant.plotter.sys_leads_hopping_pos(fsys, hop_lead_nr[0])
    N=len(hoppings[0])
    ovcurrent=[0]*2*N
    channels=len(psi)
    for wfn in psi:
        J=kwant.operator.Current(fsys)
        current=J(wfn)
        ovcurrent=list(map(operator.add, current,ovcurrent))
    values=list(ovcurrent)[:N]
    normvalues=values/max(values)
    plt.figure(figsize=(8.25,10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[19, 1]) 
    ax2 = plt.subplot(gs[0]) 
    ax1 = plt.subplot(gs[1])  
    dim=abs(xaxis[1]-xaxis[0])
    for j in range(N):
        if normvalues[j]<0:
            ax2.annotate("", xy=(hoppings[1][j][0], hoppings[1][j][1]),xytext=(hoppings[0][j][0],hoppings[0][j][1]),arrowprops=dict(color=cm.inferno(int(round(abs(values[j]*256/channels)))),headlength=100/dim,headwidth=150/dim*abs(normvalues[j]),width=70/dim*abs(normvalues[j])))
        if normvalues[j]>=0:
            ax2.annotate("", xy=(hoppings[0][j][0], hoppings[0][j][1]),xytext=(hoppings[1][j][0],hoppings[1][j][1]),arrowprops=dict(color=cm.inferno(int(round(abs(values[j]*256/channels)))),headlength=100/dim,headwidth=150/dim*abs(normvalues[j]),width=70/dim*abs(normvalues[j])))
    kwant.plotter.plot(fsys,site_size=0,site_color='white',hop_color='white',lead_color='black',lead_site_size=0.25,ax=ax2)
    ax2.set_yticks([])
    ax2.set_xticks([])
    ax2.set_aspect('equal')
    # Set the colormap and norm to correspond to the data for which the colorbar will be used.
    cmap = mpl.cm.inferno
    norm = mpl.colors.Normalize(vmin=0, vmax=channels)
    #ClorbarBase derives from ScalarMappable and puts a colorbar in a specified axes, so it has everything needed for a                                                         # standalone colorbar.  There are many more kwargs, but the following gives a basic continuous colorbar with ticks and labels
    cb1 = mpl.colorbar.ColorbarBase(ax=ax1, cmap=cmap,norm=norm,orientation='horizontal')
    cb1.set_label('Local probability current')
    plt.tight_layout()
    plt.savefig('webservice/user_static/img/current.png' ,bbox_inches='tight')


