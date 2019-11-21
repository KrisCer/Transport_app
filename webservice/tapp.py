def plot_numbers(sys, lead0, lead1, a, b, v1, v2, p1, p2, A, D, W):   #Plots the junctions with numbered scattering region atoms 
        import kwant
        from math import sqrt
        import matplotlib.pyplot as plt
        import numpy as np
        sites = list(sys.sites()) #get the sites from kwant builder, the coordinates are shown in terms of lattice vectors
        fig, ax = plt.subplots(figsize=(10,10)) 
        ax.set_aspect('equal','box')
        ax.set_yticks([])
        ax.set_xticks([])
        for count,atom in enumerate(sites): #go through all the sites in scattering region
            #Check which sublattice the atom belongs as different basis vectors p1/2 are used
            if str(atom[0])[19]=='0':
                p=p1
                xshift=0
                yshift=0.05
            else:
                p=p2
                xshift=0.05
                yshift=-0.05
            if A==120: #Checks for the angle and plots a numeber for each atom
                if atom[1][0]>-D and atom[1][0]+atom[1][1]<3: #only take the atoms in scattering centre,
                    #ignoring the lead sites
                    #atom gives the number of lattice vectors v1/2 to get to the position form (0,0)
                    #plus the basis p1/2 and a small shift for better plotting
                    ax.text(atom[1][0]*v1[0]+atom[1][1]*v2[0]+p[0]+xshift,atom[1][0]*v1[1]+atom[1][1]*v2[1]+p[1]+yshift,
                            s=str(count),color='r') 
            if A==60:
                if atom[1][0]>-sqrt(3)*1.5 and atom[1][1]<3:
                    ax.text(atom[1][0]*v1[0]+atom[1][1]*v2[0]+p[0]+xshift,atom[1][0]*v1[1]+atom[1][1]*v2[1]+p[1]+yshift,
                            s=str(count),color='r')
            if A==180:
                if atom[1][0]>=-1 and atom[1][0]<=2*D-2:
                    ax.text(atom[1][0]*v1[0]+atom[1][1]*v2[0]+p[0]+xshift,atom[1][0]*v1[1]+atom[1][1]*v2[1]+p[1]+yshift,
                            s=str(count),color='r')
        kwant.plotter.plot(sys,site_size=0,site_color='black',lead_site_size=0.25,site_edgecolor='black',lead_color='black',ax=ax) #plots the outline of the system
        xaxis=ax.get_xlim() #get the limits of x and y for further plotting
        yaxis=ax.get_ylim()
        fig.savefig('webservice/user_static/img/plot.png',bbox_inches='tight')
        return xaxis, yaxis


def make_system(A,W,F=0,D=1,t=2.75,S=0):
    import kwant
    from math import sqrt
    import matplotlib.pyplot as plt
    import numpy as np
#Check if width is even or odd and assign a variable N, that is used to construct system
    if W % 2 == 0:
        N=(W-2)/2
    else:
        N=(W-3)/2
        
################Define lattice, lead symetries, vectors and basis positions depending on angle#################

    v1,v2=(sqrt(3)/2, 0.5), (0, 1) #Lattice vectors
    if A==120:
        p1,p2=(-sqrt(3)/6,-0.5),(-sqrt(3)/3,0) #Basis positions      
        graphene = kwant.lattice.general([v1,v2],[p1,p2],norbs=1) #define graphene lattice
        a, b = graphene.sublattices #define sublatices        
        ysym = kwant.TranslationalSymmetry([sqrt(3)/2,1.5])  # lattice symmetry for lead1 in 120 degree case
        for sl in graphene.sublattices:
            ysym.add_site_family(sl, other_vectors=[(1, -1)]) # Add vector to have correct facing lead cells
        # Additinal parameter to correctly assign radius, depending on what is the distance between leads
        Z=3
        if D%2==0:
            Z=6
        r=sqrt((0.5+D/2+0.5*N)**2+(sqrt(3)/Z+(sqrt(3)/6)*N)**2)+0.01 # Asign radius to build scattering region
    #Other angles go through the same procedure, but with different vectors/radiuses
    if A==60:
        p1,p2=(sqrt(3)/3,0),(sqrt(3)/6,-0.5)       
        graphene = kwant.lattice.general([v1,v2],[p1,p2],norbs=1) 
        a, b = graphene.sublattices 
        ysym = kwant.TranslationalSymmetry([-sqrt(3)/2,1.5]) 
        r=sqrt((1.5+1*N)**2+(sqrt(3)/6)**2)+0.01
    if A==180:
        p1,p2=(sqrt(3)/3,0),(sqrt(3)/6,-0.5)
        graphene = kwant.lattice.general([v1,v2],[p1,p2],norbs=1) 
        a, b = graphene.sublattices 
        ysym = kwant.TranslationalSymmetry([-sqrt(3),0])
    if A==120 and D==2: #additional conditional to get correct looking lead connections
        xsym = kwant.TranslationalSymmetry([sqrt(3),0])
    else:
        xsym = kwant.TranslationalSymmetry([-sqrt(3),0])       
    
#######################Building system##############################

    sys = kwant.Builder() # Initialize system builder
    def circle(pos): #define the scattering region dimensions depending on angle
        x, y = pos
        if A==180:
            return -1<=x<2*sqrt(3)/3+(D-1)*sqrt(3) and -0.5<=y<=0.50+N+S #retruns the scattering region dimensions
        if A==120: #different positions depending on the distance between the leads to have correct looking lead attachments
            if D==3:
                return (x+sqrt(3)/6*(1-N)) ** 2 + (y+0.5*(1+N)) ** 2 <= r**2
            if D==2:
                return (x-sqrt(3)/6*N) ** 2 + (y+0.5*N) ** 2 <= r ** 2            
            if D==1:
                return (x-sqrt(3)/6*(1+N)) ** 2 + (y-0.5*(1-N)) ** 2 <= r**2
        if A==60:
            return (x-3*sqrt(3)/6*N) ** 2 + (y+0.5*N) ** 2 <= r ** 2 
    sys[graphene.shape(circle, (0, 0))] = 0    # Build scattering region
    sys[graphene.neighbors()] = -t #Assign nearest-neighbour hoppings
    sys.eradicate_dangling() #Delete dangling atoms
    sites=list(sys.sites()) #Obtain list of sites for further use 
    
    
    #In the case of 120 degrees for wide ribbons, the scattering region overlaps with the leads, 
    #removing atoms outside the intended scattering region
    if A==120: 
        for atom in sites:  
            if atom[1][0]<=-D: #remove extra atoms on left side
                if str(atom[0])[19]=='0':
                    del sys[a(atom[1][0],atom[1][1])] 
                else:
                    del sys[b(atom[1][0],atom[1][1])]     
            if atom[1][0]+atom[1][1]>=3: #remove extra atoms on right side
                if str(atom[0])[19]=='0':
                    del sys[a(atom[1][0],atom[1][1])] 
                else:
                    del sys[b(atom[1][0],atom[1][1])] 
    
    #Further building of the scattering center by adding lead unit cells for neater looking system and 
    #easier plotting of properties
    #Adding 2 unit cells of lead atoms to scattering region
    def lead_sites_left(pos): 
    #defining a function to build the scattering center lead cells depending on angle
    #simialr to scattering centre building previously
        x, y = pos
        if A==120:
            if W % 2 == 0:
                if F==0:
                    return -N+1.5-D/2 <= y <= 2-D/2 and -(D+5)*sqrt(3)/2<=x<=-sqrt(3)/2
                if F==1:
                    return -N+(1-D/2) <= y <= 1.5-D/2 and -(D+5)*sqrt(3)/2<=x<=-sqrt(3)/2
            else:
                return -N-(D/2-0.5) <= y <= 1.5-D/2 and -(D+5)*sqrt(3)/2<=x<=-sqrt(3)/2
        if A==60:
            if W % 2 == 0:
                if F==0:
                    return -N <= y <= 0.5 and -7*sqrt(3)/2<=x<=-sqrt(3)/2
                    
                if F==1:
                    return -N-0.5 <= y <= 0 and -7*sqrt(3)/2<=x<=-sqrt(3)/2
            else:
                return -N-0.5 <= y <= 0.5 and -7*sqrt(3)/2<=x<=-sqrt(3)/2
        if A==180:
            if W % 2 == 0:
                if F==1:
                    return 0 <= y <= 0.5+N and -5*sqrt(3)/2<=x<=-sqrt(3)/2
                    
                if F==0:
                    return -0.5 <= y <= 0+N and -5*sqrt(3)/2<=x<=-sqrt(3)/2
            else:
                return -0.5 <= y <= 0.5+N and -5*sqrt(3)/2<=x<=-sqrt(3)/2
    sys[graphene.shape(lead_sites_left, (-sqrt(3)/2*D, 0))] = 0    #adds the scattering centre lead sites
    sys[graphene.neighbors()] = -t #adds hoppings between sites
    
    #Adding lead cells to the other side, which is more tideous as the leads go in an angle and are not easy to specify
    #hence we build a cell by indivudually placing atoms and then repeat the cell both lengthwise and widthwise
    if A==120:
        for Z in range(3): #adds 3 unit cells
            sys[a(1+1*Z,1+1*Z)]=0 #manually adding an atom at specified coordinate with onsite potential of 0
            sys[a(1+1*Z,2+1*Z)]=0
            sys[a(2+1*Z,1+1*Z)]=0
            sys[b(1+1*Z,1+1*Z)]=0    
            sys[b(2+1*Z,0+1*Z)]=0
            sys[b(2+1*Z,1+1*Z)]=0
            if N>0: #Checking if widht>3
                N=int(N)
                for K in range(N): #Extra atoms added to account for width as original builds 3-AGNR leads
                    K=K+1
                    sys[a(1+K+1*Z,1-K+1*Z)]=0 #depending on the width parameter N, extend the lead widths
                    sys[a(2+K+1*Z,1-K+1*Z)]=0 
                    sys[b(2+K+1*Z,1-K+1*Z)]=0
                    sys[b(2+K+1*Z,0-K+1*Z)]=0
            if N==0:
                K=0
        if W % 2 == 0: #Checking if lead width is even, hence posibility of 2 chiralities 
            #Need to compleatly rework the atom positions for one of the chiralities by manually adding/removing 
            #atoms from the initial "cell"
            if F==0:
                del sys[a(2+K,1-K)] #deleting a specific atom 
                del sys[a(3+K,2-K)]
                del sys[b(3+K,1-K)]
                del sys[a(4+K,3-K)]
                del sys[b(4+K,2-K)]
                del sys[a(2+K,2-K)]
                del sys[a(3+K,3-K)]
                del sys[b(2+K,1-K)]
                del sys[b(3+K,2-K)]
                del sys[b(4+K,3-K)]
                sys[a(1,3)]=0
                sys[a(2,4)]=0
                sys[b(1,2)]=0
                sys[a(3,4)]=0
                sys[b(2,3)]=0
                sys[b(3,4)]=0             
            if F==1:
                del sys[a(2+K,1-K)]
                del sys[a(3+K,2-K)]
                del sys[b(3+K,1-K)]
                del sys[a(4+K,3-K)]
                del sys[b(4+K,2-K)]
    if A==60: #similar to 120 case
        for Z in range(3):
            sys[a(-2-1*(Z-1),4+2*(Z-1))]=0
            sys[a(-2-1*Z,3+2*Z)]=0
            sys[a(-3-1*(Z-1),4+2*(Z-1))]=0
            sys[b(-2-1*(Z-1),4+2*(Z-1))]=0    
            sys[b(-2-1*Z,3+2*Z)]=0
            sys[b(-1-1*Z,3+2*Z)]=0
            if N>0:
                N=int(N)
                for K in range(N):
                    K=K+1
                    sys[a(-2+K-1*Z,3+2*Z)]=0
                    sys[a(-2+K-1*(Z-1),4+2*(Z-1))]=0
                    sys[b(-1+K-1*Z,3+2*Z)]=0    
                    sys[b(-2+K-1*(Z-1),4+2*(Z-1))]=0
            if N==0:
                K=0
            if W % 2 == 0: 
                #more streamlined process than the 120 case with less lines
                #to account for chiraliry of even width leads
                if F==0:
                    if Z!=0:
                        del sys[a(-2-(Z-1)+K,4+2*(Z-1))]
                    del sys[b(-1-Z+K,3+2*Z)]
                if F==1:
                    if Z!=0:
                        del sys[a(-3-(Z-1),4+2*(Z-1))]
                    del sys[b(-2-Z,3+2*Z)]
                    
    def lead_sites_right(pos): 
        #180 case allows us to build a simple function once again to add the leads sites to the
        #other side of the scattering region. simply repeating the process from lead_sites_left but adding the shift S
        x, y = pos
        if W % 2 == 0:
            if F==0:
                return +S <= y <= 0.5+S+N and (D-0.5)*sqrt(3)<=x<=(D+2-0.5)*sqrt(3)
            if F==1:
                return -0.5+S <= y <= 0+S+N and (D-0.5)*sqrt(3)<=x<=(D+2-0.5)*sqrt(3)
        else:
            return -0.5+S <= y <= 0.5+S+N and (D-0.5)*sqrt(3)<=x<=(D+2-0.5)*sqrt(3)
        
    if A==180:
        sys[graphene.shape(lead_sites_right, (D*sqrt(3),0.5+S))] = 0
        
    sys[graphene.neighbors()] = -t
    
############################### Lead building ###################################################
    #Building the unit cell of the left lead, similar process to building the lead cells in
    #scattering centre, but this time without specifying x component as it will be semi-infinite
    #and will have a symetry 
    def lead0_shape(pos): 
        x, y = pos
        if A==120:
            if W % 2 == 0:
                if F==0:
                    return -N+1.5-D/2 <= y <= 2-D/2
                if F==1:
                    return -N+(1-D/2) <= y <= 1.5-D/2
            else:
                return -N-(D/2-0.5) <= y <= 1.5-D/2 
        if A==60:
            if W % 2 == 0:
                if F==0:
                    return -N <= y <= 0.5
                    
                if F==1:
                    return -N-0.5 <= y <= 0
            else:
                return -N-0.5 <= y <= 0.5
            
        if A==180:
            if W % 2 == 0:
                if F==1:
                    return 0 <= y <= 0.5+N
                    
                if F==0:
                    return -0.5 <= y <= 0+N
            else:
                return -0.5 <= y <= 0.5+N

    lead0 = kwant.Builder(xsym) #starting the lead builder and specifying the symetry of semi-infinite direction
    lead0[graphene.shape(lead0_shape, (0,0.25))] = 0
    lead0[graphene.neighbors()] = -t
       
    
    #Building the other lead at some angle. Specifying the lead cell by manually adding atoms, that
    #have translational symetry ysym and will be basis of semi-infinite lead.
    lead1 = kwant.Builder(ysym)
    def lead1_shape(pos):
        x, y = pos
        if W % 2 == 0:
            if F==0:
                return +S <= y <= 0.5+S+N
            if F==1:
                return -0.5+S <= y <= 0+S+N
        else:
            return -0.5+S <= y <= 0.5+S+N
    if A==180:
        lead1[graphene.shape(lead1_shape, (0,0.5+S))] = 0
    if A==60:
        lead1[a(-2,4)]=0
        lead1[a(-2,3)]=0
        lead1[a(-3,4)]=0
        lead1[b(-2,4)]=0    
        lead1[b(-2,3)]=0
        lead1[b(-1,3)]=0
        if N>0:
            N=int(N)
            for K in range(N):
                K=K+1
                lead1[a(-2+K,3)]=0
                lead1[a(-2+K,4)]=0
                lead1[b(-1+K,3)]=0    
                lead1[b(-2+K,4)]=0
        if N==0:
            K=0
        if W % 2 == 0:
            if F==0:
                del lead1[a(-2+K,4)]
                del lead1[b(-1+K,3)]
            if F==1:
                del lead1[a(-3,4)]
                del lead1[b(-2,3)]
    if A==120:
        lead1[a(1,2)]=0
        lead1[a(2,1)]=0
        lead1[a(2,2)]=0
        lead1[b(2,2)]=0    
        lead1[b(2,1)]=0
        lead1[b(3,1)]=0
        if N>0:
            N=int(N)
            for K in range(N):
                K=K+1
                lead1[a(2+K,1-K)]=0
                lead1[a(2+K,2-K)]=0  
                lead1[b(2+K,1-K)]=0
                lead1[b(3+K,1-K)]=0
        if N==0:
            K=0
        if W % 2 == 0:
            if F==0:
                del lead1[a(2+K,1-K)]
                del lead1[a(1+K,1-K)]
                del lead1[b(2+K,0-K)] 
                del lead1[b(2+K,1-K)]
                lead1[a(0,2)]=0
                lead1[b(1,2)]=0 
            if F==1:
                del lead1[a(2+K,1-K)]
                del lead1[b(2+K,0-K)]
    
    lead1[graphene.neighbors()] = -t
    
    if A==120 and D==2: #depending on the angle/width attach the first lead to obtain correct display
        sys.attach_lead(lead0.reversed()) 
    else:    
        sys.attach_lead(lead0)
    if A==180: #depending on the angle attach 2nd lead for correct transitions
        sys.attach_lead(lead1.reversed())
    else:
        sys.attach_lead(lead1)
     
    #assigns global variables that are used for further calculations    
    return sys, lead0, lead1, a, b, v1, v2, p1, p2, A, D, W
