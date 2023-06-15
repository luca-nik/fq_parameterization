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
