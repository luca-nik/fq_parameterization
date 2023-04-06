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
        self.number_connections = [] 
        self.connected_to = [] 
        self.surface_atoms = [] 

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
    def get_connectivity(self, print_info = False):
        #
        """Procedure to get the connectivity information of each atom of the molecule"""
        #
        # sanity check
        #
        if (self.atoms == 0):
            print('ERROR: you need to initalize the molecule before computing the connectivity of its atoms')
            sys.exit()
        #
        threshold = 1.5 #in Angstrom
        #
        connection_information = []
        #
        if print_info:
            print('Printing the connectivity information:')
        #
        # Get the connetivity information of each atom according to the threshold distance
        #
        for i,at in enumerate(self.atomtypes):
            connected_to = []
            connections = 0
            for j,at2 in enumerate(self.atomtypes):
                if (i != j):
                    dist = np.linalg.norm(self.coords[i] - self.coords[j])
                    if (dist <= threshold):
                       connected_to.append(j)
                       connections += 1
            #
            self.number_connections.append(connections)
            self.connected_to.append(connected_to)
            #
            # Print information abut the connectivity if required
            #
            if print_info:
                info_string = ''
                #
                for index, connector in enumerate(self.connected_to[i]):
                    info_string += str(self.atomtypes[connector]) + ' ' + str(connector)
                    if index < len(self.connected_to[i])-1:
                        info_string += ', '
                #
                print(at + ' ' + str(i).ljust(len(str(self.atoms))) +  ' - connected to ' + str(self.number_connections[i]).ljust(2) + ' atoms : ' + info_string)
    #
    def get_interface_atoms(self,print_info = False):
        #
        """Procedure to get the atoms which are at the interface with the solvent""" ## DA SPIEGARE MEGLIO STA COSA
        #
        # Sanity check
        #
        if (len(self.number_connections) == 0):
            print('ERROR: you need to compute the connectivities before identifying interface atoms')
            sys.exit()
        #
        if print_info:
            print('Printing the atoms which are at the interface with the solvent:')
        #
        number_connections_expected =  constants.number_connections()
        #
        for i, at in enumerate(self.atomtypes):
            if self.number_connections[i] <= number_connections_expected[at]:
                self.surface_atoms.append(i)
            #
            if print_info:
               print(at + ' ' + str(i).ljust(len(str(self.atoms)))) 
    #

