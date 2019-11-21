#Used to remove specified atoms form the scattering centre by entering a list of atom numbers
def del_atoms(sys,a,b,atoms):
    import kwant
    atoms=sorted(atoms, reverse=True) #sort so it starts from the largest one, such that the sites output doesn't change
    sites = list(sys.sites()) #retrieve the list of atoms and their positions in terms of the lattice vecotrs
    for atom in atoms:
        vec1=int(sites[atom][1][0]) #get the number of 1st vector
        vec2=int(sites[atom][1][1]) #numebr of 2nd vector
        lattice=str(sites[atom][0])[19] #obtain the sublattice 
        if lattice=='0':
            del sys[a(vec1,vec2)] #delete the atom
        else:
            del sys[b(vec1,vec2)]            
    return sys #return the updated system
