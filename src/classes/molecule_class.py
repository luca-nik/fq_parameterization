import numpy as np
import constants
import sys
#
class molecule:
    #
    """Molecule class object"""
    #
    def __init__(self, atomtypes=[],coordinates=[]):
        #
        """Initialization of the molecule class object"""
        #
        self.atomtypes = list(atomtypes).copy()
        self.atoms = len(self.atomtypes)
        self.coords = np.asarray(coordinates.copy())
        self.molecules = [0 for i in range(self.atoms)]

    def write_xyz(self, name = 'nome', directory = '', comment = ''):
        #
        """Procedure to write the molecule object into a formatted xyz file"""
        #
        if directory[-1] != '/':
            directory +='/'
        with open(directory+name+'.xyz', 'w') as outfile:
            outfile.write(str(self.atoms) + '\n')
            outfile.write(str(comment))
            for i,sym in enumerate(self.atomtypes):
                outfile.write('\n' + sym.rjust(2) + '  ' + \
                               '{:5.5f}'.format(self.coords[i][0]).rjust(10)+ '  ' + \
                               '{:5.5f}'.format(self.coords[i][1]).rjust(10)+ '  ' + \
                               '{:5.5f}'.format(self.coords[i][2]).rjust(10))
    #
    def initialize_from_xyz(self, file):
        #
        """Procedure to initialize a molecule class object from a xyz file"""
        #
        with open(file, 'r') as f:
            lines = f.readlines()
        self.atoms = int(lines[0].split()[0])
        coords = []
        attypes = []
        for line in lines[2:]: 
            if len(line.split()) != 0:
                coords.append([float(i) for i in line.split()[1:]])
                self.atomtypes.append(line.split()[0])
        self.coords = np.asarray(coords)
        self.molecules = [0 for i in range(self.atoms)]
    #
    def get_bond_axes(self):
        #
        """Procedure to get the axes along which you have the bonds"""
        #
        for i,at in enumerate(self.atomtypes):
            pass
            #constants.
