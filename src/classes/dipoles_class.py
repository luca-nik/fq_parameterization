import numpy as np
import constants
import sys
#
class dipoles:
    #
    """Dipoles class object"""
    #
    def __init__(self, n_dipoles = 0, positions = [], directions= []):
        #
        """Initialization of the dipoles class object"""
        """n_dipoles  = number of dipoles;
           positions  = position of the dipole center;
           directions = versor of the direction towards which the dipole was created""" 
        #
        self.n_dipoles  = n_dipoles
        self.positions  = np.asarray(positions.copy())
        self.directions = np.asarray(directions.copy())

    def write_xyz(self, name = 'nome', directory = './', comment = ''):
        #
        """Procedure to write the dipole object into a formatted xyz file"""
        #
        if directory[-1] != '/':
            directory +='/'
        #
        with open(directory+name+'.xyz', 'w') as outfile:
            outfile.write(str(self.n_dipoles) + '\n')
            outfile.write(str(comment))
            for i in range(0,self.n_dipoles):
                outfile.write('\n' + 'Mg' + '  ' + \
                               '{:5.5f}'.format(self.positions[i][0]).rjust(10)+ '  ' + \
                               '{:5.5f}'.format(self.positions[i][1]).rjust(10)+ '  ' + \
                               '{:5.5f}'.format(self.positions[i][2]).rjust(10))
    #
    def write_dip(self, name = 'nome', directory = './'):
        #
        """Procedure to write the .dip file needed to move the dipoles in the space"""
        #
        if directory[-1] != '/':
            directory +='/'
        #
        with open(directory+name+'.dip', 'w') as outfile:
            outfile.write('total number of dipoles: ' + str(self.n_dipoles) + '\n')
            outfile.write('Dipole | ' + 'Pos X'.rjust(10) + ' ' +'Pos Y'.rjust(10) + ' ' + 'Pos Z'.rjust(10) + \
                           ' | ' + 'Dir X'.rjust(10) + ' ' +'Dir Y'.rjust(10) + ' ' + 'Dir Z'.rjust(10))
            for i in range(0,self.n_dipoles):
                outfile.write('\n' + str(i).ljust(6) + ' | ' + \
                               '{:5.5f}'.format(self.positions[i][0]).rjust(10)+ ' '   + \
                               '{:5.5f}'.format(self.positions[i][1]).rjust(10)+ ' '   + \
                               '{:5.5f}'.format(self.positions[i][2]).rjust(10) + ' | '+ \
                               '{:5.5f}'.format(self.directions[i][0]).rjust(10)+ ' '  + \
                               '{:5.5f}'.format(self.directions[i][1]).rjust(10)+ ' '  + \
                               '{:5.5f}'.format(self.directions[i][2]).rjust(10))
    #
    def initialize_from_dip(self, file):  
        #                                  
        """Procedure to initialize a dipoles object from a .dip file"""
        #
        with open(file, 'r') as f:
            lines = f.readlines()
        self.n_dipoles = int(lines[0].split(':')[1])
        positions = []
        directions = []
        for line in lines[2:]: 
            if len(line.split()) != 0:
                coordinates = line.split('|')[1]
                positions.append([float(i) for i in coordinates.split()[:]])
                versors = line.split('|')[2]
                directions.append([float(i) for i in versors.split()[:]])
        #
        self.positions  = np.asarray(positions)
        self.directions = np.asarray(directions)
    #
#    def move_dipoles(self, which_dipoles = [], displacements = []):
#        #
#        """Procedure to move some dipoles of a certain displacement"""
#        """which_dipoles is a list of integers identyfing the dipoles
#           displacements is a list of floats which are expressed in Angstrom
#           the displacement happen along the line and the direction identified by dipole.directions
#        """
#        #

                                          
                                          
                                          
