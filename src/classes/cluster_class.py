import numpy as np
import constants
import sys
from classes import molecule_class
#
class cluster:
    #
    """Molecule class object"""
    #
    def __init__(self, molecules = [], mol_id = 'OHH'):
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
        #
        if (type(mol_id) != str):
            print('ERROR: initializing a cluster object with wrong mol_id to identify the molecules')
            sys.exit()
        #
        if (mol_id == ''):
            print('ERROR: initializing a cluster object with empty mol_id to identify the molecules')
            sys.exit()
        #
        self.mol_id = mol_id


    def write_clust(self, name = 'nome', directory = './'):
        #
        """Procedure to write the cluster object into a formatted clust file"""
        #
        if directory[-1] != '/':
            directory +='/'
        with open(directory+name+'.clust', 'w') as outfile:
            #
            outfile.write('Molecules identifier: ' + self.mol_id + ' \n')
            #
            for mol_indx, molecule in self.molecules:
                outfile.write('molecule :' + str(mol_indx) + ', atoms: ' + str(molecule.atoms) + '\n')
                if (molecule.charge != 0):
                    outfile.write('  CHARGE:' + str(charge))
                for i,sym in enumerate(molecule.atomtypes):
                    outfile.write('\n' + sym.rjust(2) + '  ' + \
                                   '{:5.5f}'.format(molecule.coords[i][0]).rjust(10)+ '  ' + \
                                   '{:5.5f}'.format(molecule.coords[i][1]).rjust(10)+ '  ' + \
                                   '{:5.5f}'.format(molecule.coords[i][2]).rjust(10))
    #
    #
    #
    #def initialize_from_clust(self, file):
    # #   #
    #    """Procedure to initialize a cluster class object from a clust file"""
    #    #
    #    with open(file, 'r') as f:
    #        lines = f.readlines()
    #    #

    #    self.atoms = int(lines[0].split()[0])
    #    coords = []
    #    attypes = []
    #    #
    #    # Get the charge
    #    #
    #    try:
    #        self.charge = float(lines[1].split('CHARGE:')[-1])
    #    except:
    #        pass
    #    #
    #    # Get the coordinates
    #    #
    #    for line in lines[2:]: 
    #        if len(line.split()) != 0:
    #            coords.append([float(i) for i in line.split()[1:]])
    #            self.atomtypes.append(line.split()[0])
    #    self.coords = np.asarray(coords)
    #    self.molecules = [0 for i in range(self.atoms)]

