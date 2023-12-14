import numpy as np
import constants
import sys
from classes import molecule_class
#
class cluster:
    #
    """Molecule class object"""
    #
    def __init__(self, molecules = []):
        #
        """Initialization of the cluster class object"""
        #
        if molecules != []:
            for molecule in molecules:
                if (not isinstance(molecule, molecule_class.molecule)):
                    print('ERROR: initializing a cluster object with non-molecule elements')
                    sys.exit()
        #
        self.molecules = molecules.copy()
        ##
        #if (type(mol_id) != str):
        #    print('ERROR: initializing a cluster object with wrong mol_id to identify the molecules')
        #    sys.exit()
        ##
        #if (mol_id == ''):
        #    print('ERROR: initializing a cluster object with empty mol_id to identify the molecules')
        #    sys.exit()
        ##
        #self.mol_id = mol_id


    def write_clust(self, name = 'nome', directory = './'):
        #
        """Procedure to write the cluster object into a formatted clust file"""
        #
        if directory[-1] != '/':
            directory +='/'
        with open(directory+name+'.clust', 'w') as outfile:
            #
            for mol_indx, molecule in enumerate(self.molecules):
                #
                if (molecule.charge != 0):
                    outfile.write('molecule: ' + str(mol_indx+1).ljust(5) + '; atoms: ' + str(molecule.atoms).ljust(5) + \
                                  '; charge: ' + str(molecule.charge).ljust(5) + '\n')
                else:
                    outfile.write('molecule: ' + str(mol_indx+1).ljust(5) + '; atoms: ' + str(molecule.atoms).rjust(5) + '\n')
                   
                #
                for i,sym in enumerate(molecule.atomtypes):
                    outfile.write(sym.rjust(2) + '  ' + \
                                  '{:5.5f}'.format(molecule.coords[i][0]).rjust(10)+ '  ' + \
                                  '{:5.5f}'.format(molecule.coords[i][1]).rjust(10)+ '  ' + \
                                  '{:5.5f}'.format(molecule.coords[i][2]).rjust(10) + '\n')
                outfile.write('\n')
    #
    #
    #
    def initialize_from_clust(self, file):
        #
        """Procedure to initialize a cluster class object from a clust file"""
        #
        with open(file, 'r') as f:
            lines = f.readlines()
        #
        cluster_molecules = []
        indices = [i for i,line in enumerate(lines) if 'molecule:' in line]
        #
        for index in indices:
            #
            molecule = molecule_class.molecule()
            atoms = lines[index].split(';')[1]
            molecule.atoms = int(atoms.split(':')[-1])
            #
            if ('charge:' in lines[index]):
                charge = lines[index].split(';')[-1]
                molecule.charge = float(charge.split(':')[-1])
            #
            coords = []
            at_types = []
            #
            for i,line in enumerate(lines[index+1:index+1+molecule.atoms]):
                c = [float(i) for i in line.split()[1:]]
                coords.append(c)
                at_t = line.split()[0]
                at_types.append(at_t)
            #
            molecule.coords = coords
            molecule.atomtypes = at_types
            cluster_molecules.append(molecule)
        #
        self.molecules = cluster_molecules.copy()
            
    #
    def get_atomtypes(self):
        #
        """Procedure to get a list of atomtypes of the cluster object"""
        atomtypes = []
        for mol in self.molecules:
            for i in mol.atomtypes:
                if (i not in atomtypes):
                    atomtypes.append(i)
        #
        return atomtypes
    #
    def get_rdf(self, molecule_1 = '', molecule_2 = '', r_max = 30, bins = 1000):
        #
        # r max in angstorm
        #
        if isinstance(molecule_1, int) and isinstance(molecule_2, str):
            if molecule_2 == '':
                print('RDF but molecule_2 not set')
                sys.exit()
            else:
                dists = []
                target_index = molecule_1 - 1 
                target = self.molecules[target_index]
                #
                coordinates = target.get_cm()
                coordinates = np.array(coordinates)
                
                for index, molecule in enumerate(self.molecules):
                    if index == target_index:
                        pass
                    else: #If a system has more atoms of molecule_2 atomtypes, then you take only the one with smaller distance
                        at_indices = [at_index for at_index, atom in enumerate(molecule.atomtypes) if atom == molecule_2]
                        mindist = -1
                        for atom_index in at_indices:
                            coordinates2 = np.array(molecule.coords[atom_index])
                            #
                            dist = np.linalg.norm(coordinates-coordinates2)
                            if atom_index == at_indices[0]:
                                mindist = dist
                            else:
                                if mindist > dist:
                                    mindist = dist
                        #
                        if mindist <= r_max:
                            dists.append(mindist)

                #for index, molecule in enumerate(self.molecules):
                #    if index == target_index:
                #        pass
                #    else: #only first occurrence
                #        atom_index = molecule.atomtypes.index(molecule_2)
                #        coordinates2 = np.array(molecule.coords[atom_index])
                #        #
                #        for coordinates in target.coords:
                #            #
                #            dist = np.linalg.norm(np.array(coordinates)-coordinates2)
                #            if dist <= r_max:
                #                dists.append(dist)
        #
        elif isinstance(molecule_1, str) and isinstance(molecule_2, str):
            if molecule_2 == '' or molecule_1 == '':
                print('RDF but molecule_2 or molecule_1 not set')
                sys.exit()
            else:
                dists = []
                #
                for index, molecule2 in enumerate(self.molecules):
                    atom_index2 = molecule2.atomtypes.index(molecule_2)
                    coordinates2 = np.array(molecule2.coords[atom_index2])
                    #
                    for iindex, molecule1 in enumerate(self.molecules):
                        if iindex == index:
                            pass
                        else:
                            atom_index1 = molecule1.atomtypes.index(molecule_1)
                            coordinates1 = np.array(molecule1.coords[atom_index1])

                            dist = np.linalg.norm(coordinates1-coordinates2)
                            dists.append(dist)


        else:
            raise ValueError("Input must be an integer or a string.")
        #
        dists = np.array(dists)

        # Create bin edges
        bin_edges = np.linspace(0, r_max, bins + 1)

        # Populate bins
        hist, edges = np.histogram(dists, bins=bin_edges, density = False)
        #
        # Normalize RDF
        #
        system_volume = (4/3)*np.pi*np.amax(dists)**3
        #
        N_particles = (len(dists))#/target.atoms) #valido solo se prendo un atomo che c'e' solo una volta nelle molecole
        #
        ##Number density
        rho = N_particles / system_volume
        #print(rho)

        ## Calculate the volume of each bin
        r_inner = edges[:-1]
        r_outer = edges[1:]
        #
        shell_volume = (4/3) * np.pi * (r_outer**3 - r_inner**3)
        #
        ## Normalize the RDF
        bin_midpoints = (edges[1:] + edges[:-1]) / 2
        rdf = hist / (4*np.pi*rho*edges[1:]**2*r_max / bins)#*target.atoms)
        test = rdf*bin_midpoints**2
        rdf_integral = np.sum(test*r_max/bins)*4*np.pi*rho#np.rect(test, bin_midpoints)*4*np.pi*rho

        #print("Integral of the RDF:", rdf_integral)
        #sys.exit()
        #
        return rdf, bin_midpoints

