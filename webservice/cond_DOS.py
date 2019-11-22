#Plots the conductance from emin to emax (in terms of t) in 300 steps.
def plot_conductance(fsys,perfsys,emin=-3, emax=3, t=2.75): #define what to calculate for the structures
    import numpy as np
    import kwant
    import matplotlib.pyplot as plt
    trans=[] #empty list to store values
    trans2=[]
    energies=np.linspace((emin+0.0001)*t, (emax+0.0001)*t, num=100) #specify energy grid
    for energy in energies:
        smatrix = kwant.smatrix(fsys, energy,check_hermiticity=False) #calculates the S matrix of the system in interest
        trans.append(smatrix.transmission(1,0)) #calcualtes the transmission coefficient
        smatrix2 = kwant.smatrix(perfsys, energy,check_hermiticity=False) #same but for ideal lead system
        trans2.append(smatrix2.transmission(1,0))
    for j in range(len(energies)):
        energies[j]=energies[j]/t #convert energies back in terms of t                        
            #Plotting comands           
    plt.figure(figsize=(10,10))
    plt.xlabel("Energy (t)")
    plt.ylabel("Conductance (e\N{SUPERSCRIPT TWO}/h)")
    plt.xticks(np.arange(emin, emax+(emax-emin)/10, round((emax-emin)/10,2)),rotation=90)
    plt.grid(b=None, which='major', axis='both')
    plt.plot(energies, trans,label='Junction')
    plt.plot(energies, trans2,label='GNR')
    plt.legend()
    plt.savefig('webservice/user_static/img/cond.png',bbox_inches='tight')

#builds a perfect GNR system to compare with the system in interest
def perfect_system(W,t=2.75): #setting up the system with vectors, lattice and symetries
    import kwant
    from math import sqrt
    p1,p2=(sqrt(3)/3,0),(sqrt(3)/6,-0.5)
    v1,v2=(sqrt(3)/2, 0.5), (0, 1)
    graphene = kwant.lattice.general([v1,v2],[p1,p2],norbs=1) 
    a, b = graphene.sublattices 
    xsym = kwant.TranslationalSymmetry([-sqrt(3),0])
    leads = kwant.Builder() # Initialize system builder
    def rec(pos):
        x, y = pos
        return -sqrt(3)<=x<=3*sqrt(3) and 0<=y<=(W-1)/2
    
    leads[graphene.shape(rec, (0, 0))] = 0    # Build scattering region
    leads[graphene.neighbors()] = -t #Assign nearest-neighbour hoppings
    leads[a(0, round((W-1)/2+1, 1))] = 0 #add two "fake" atoms for color map normalisation
    leads[b(0, -1)] = 0
    def lead_shape(pos):  #desing the leads, which are the same as scattering region, but without the x restirction 
        x, y = pos
        return 0 <= y <= (W-1)/2
    lead = kwant.Builder(kwant.TranslationalSymmetry([sqrt(3),0]))
    lead[graphene.shape(lead_shape, (0,0))] = 0
    lead[graphene.neighbors()] = -t                                                                                         
    leads.attach_lead(lead.reversed()) #attach the leads
    leads.attach_lead(lead)
    leads = leads.finalized() #finalize the system
    return leads

#Plots the DOS from emin to emax (in terms of t) in 300 steps. Same procedure as get_conductance 
def get_DOS(fsys,perfsys,emin=-3, emax=3, t=2.75):
    import numpy as np
    import kwant
    import matplotlib.pyplot as plt
    dos=[]
    dos2=[]
    energies=np.linspace((emin+0.0001)*t, (emax+0.0001)*t, num=100)
    for J in energies:
        dos.append(np.sum(kwant.ldos(fsys,J))) #sum over all atom ldos 
        dos2.append(np.sum(kwant.ldos(perfsys,J)))
    for j in range(len(energies)):
        energies[j]=energies[j]/t        
    plt.figure(figsize=(10,10))
    plt.plot(energies, dos,label='Junction')
    plt.plot(energies, dos2,label='GNR')
    plt.xlabel("Energy (t)")
    plt.ylabel("DOS")
    plt.xticks(np.arange(emin, emax+(emax-emin)/10, round((emax-emin)/10,2)),rotation=90)
    plt.grid(b=None, which='major', axis='both')
    plt.legend()
    plt.savefig('webservice/user_static/img/DOS.png',bbox_inches='tight')
        




