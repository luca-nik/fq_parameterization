import numpy as np
import constants
import sys
#
class dipoles:
    #
    """Dipoles class object"""
    #
    def __init__(self, n_dipoles = 0, positions = [], directions= [], signs = []):
        #
        """Initialization of the dipoles class object"""
        """n_dipoles  = number of dipoles;
           positions  = position of the dipole center;
           directions = versor of the direction towards which the dipole was created""" 
        #
        self.n_dipoles  = n_dipoles
        self.positions  = np.asarray(positions.copy())
        self.directions = np.asarray(directions.copy())
        #
        # Sanity checks and assignement on the signs option
        #
        self.check_and_assign_signs(signs)
        self.name = ''
    #
    # Write xyz procedure
    #
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
    # Write dip procedure
    #
    def write_dip(self, name = 'nome', directory = './'):
        #
        """Procedure to write the .dip file needed to move the dipoles in the space"""
        #
        # Sanity check
        #
        if directory[-1] != '/':
            directory +='/'
        #
        with open(directory+name+'.dip', 'w') as outfile:
            outfile.write('total number of dipoles: ' + str(self.n_dipoles) + '\n')
            outfile.write('Dipole | ' + 'Pos X'.rjust(10) + ' ' +'Pos Y'.rjust(10) + ' ' + 'Pos Z'.rjust(10) + \
                           ' | ' + 'Dir X'.rjust(10) + ' ' +'Dir Y'.rjust(10) + ' ' + 'Dir Z'.rjust(10)      + \
                           ' | ' + 'Sign'.rjust(6))
            for i in range(0,self.n_dipoles):
                outfile.write('\n' + str(i).ljust(6) + ' | ' + \
                               '{:5.5f}'.format(self.positions[i][0]).rjust(10)+ ' '    + \
                               '{:5.5f}'.format(self.positions[i][1]).rjust(10)+ ' '    + \
                               '{:5.5f}'.format(self.positions[i][2]).rjust(10) + ' | ' + \
                               '{:5.5f}'.format(self.directions[i][0]).rjust(10)+ ' '   + \
                               '{:5.5f}'.format(self.directions[i][1]).rjust(10)+ ' '   + \
                               '{:5.5f}'.format(self.directions[i][2]).rjust(10)+ ' | ' + \
                               self.signs[i].rjust(6)
                               )
    #
    # Initialize from dip procedure
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
        signs = []
        #
        # Read the information and store them
        #
        for line in lines[2:]: 
            if len(line.split()) != 0:
                coordinates = line.split('|')[1]
                positions.append([float(i) for i in coordinates.split()[:]])
                versors = line.split('|')[2]
                directions.append([float(i) for i in versors.split()[:]])
                signs.append(line.split('|')[3])
        #
        self.positions  = np.asarray(positions)
        self.directions = np.asarray(directions)
        self.check_and_assign_signs(signs) 
        self.name = file
    #
    def move_dipoles(self, which_dipoles = [], displacements = [], create_new_dipoles = True):
        #
        """Procedure to move some dipoles of a certain displacement

           which_dipoles is a list of integers identyfing the dipoles
           displacements is a list of floats which are expressed in Angstrom
           the displacement happen along the line and the direction identified by dipole.directions

           REMEMBER: the dipoles are generated already with a displacement of 1 Angstrom from
                     the reference atoms. This will add to your displacement
        """
        #
        # Sanity checks
        #
        if type(which_dipoles) != list:
            which_dipoles = [which_dipoles]
        #
        if (len(which_dipoles) == 0):
            which_dipoles = [i for i in range(0, self.n_dipoles)] #displace them all
        #
        elif(len(which_dipoles) != len(set(which_dipoles))):
            print('ERROR: move_dipoles you have duplicates in which_dipoles' )
            sys.exit()
        #
        for i in which_dipoles:
            if type(i) != int:
                print('ERROR: move_dipoles called with non-integer which_dipoles list')
                sys.exit()
        #
        if type(displacements) != list:
            displacements = [displacements]
        #
        if (len(displacements) == 0):
            print('WARNING: move_dipoles called without any displacement, dipoles object was not modified')
            return
        try:
            checked_displ = [float(i) for i in displacements]
        except ValueError:
            print('ERROR: move_dipoles called with non-float displacements')
            sys.exit()
        #
        if (len(displacements) == 1):
            checked_displ = [float(displacements[0]) for i in range(0,len(which_dipoles))] #you move all dipoles of the same amount
        else:
            if(len(displacements) != len(which_dipoles)):
               print('ERROR: move_dipoles displacements do not match the number of dipoles')
               sys.exit()
        #
        for i in which_dipoles:
            if type(i) != int:
                print('ERROR: move_dipoles called with non-integer which_dipoles list')
                sys.exit()
            elif i > self.n_dipoles:
                print('ERROR: move_dipoles moving non-existing dipole number ' + str(i)) 
                sys.exit()
        #
        # Move the dipoles
        #
        new_positions  = self.positions.copy()
        #
        for i,dip in enumerate(which_dipoles):
            new_positions[dip,:] += checked_displ[i]*self.directions[dip,:]
        #
        if create_new_dipoles:
            new_directions = self.directions.copy()
            new_dipoles = dipoles(n_dipoles = self.n_dipoles,\
                          positions = new_positions, directions = new_directions)
            return new_dipoles
        else:
            self.positions = new_positions.copy()
    #
    # change signs procedure
    #
    def change_sign(self, which_dipoles = [], signs = []):
        #
        """Procedure to change the sign of which_dipoles according to the input signs list"""
        #
        # Sanity checks
        #
        if (type(which_dipoles) != list):
            which_dipoles = [which_dipoles]
        #
        if (type(signs) != list):
            signs = [signs]
        #
        if ((len(which_dipoles) == 0 and len(signs) > 1) or \
            (len(which_dipoles) != 0 and len(signs) > len(which_dipoles))):
            print('ERROR: change_signs called but len(signs) > len(which_dipoles)')
            sys.exit()
        #
        elif (len(which_dipoles) == 0 and len(signs) == 0):
            return
        #
        elif (len(which_dipoles) == 0 and len(signs) == 1):
            signs_to_check = [signs[0] for i in range(0,self.n_dipoles)]
            self.check_and_assign_signs(signs_to_check)
        #
        elif(len(which_dipoles) != len(set(which_dipoles))):
            print('ERROR: change_signs you have duplicates in which_dipoles' )
            sys.exit()
        #
        elif (len(which_dipoles) != 0 and len(which_dipoles) == len(signs)):
            #
            for i in which_dipoles:
                if type(i) != int:
                    print('ERROR: change_signs called with non-integer which_dipoles list')
                    sys.exit()
                elif i > self.n_dipoles:
                    print('ERROR: change_signs changing sign to non-existing dipole number ' + str(i)) 
                    sys.exit()
            #
            signs_to_check = self.signs.copy()
            for index, dipole in enumerate(which_dipoles):
                signs_to_check[dipole] = signs[index]
            self.check_and_assign_signs(signs_to_check)
        #
        elif (len(which_dipoles) != 0 and len(which_dipoles) != len(signs)):
            #
            if (len(signs) == 0 or len(signs) > 1 ):
                print('ERROR: change_signs called with unclear sign assignement')
                sys.exit()
            #
            else:
                #
                for i in which_dipoles:
                    if type(i) != int:
                        print('ERROR: change_signs called with non-integer which_dipoles list')
                        sys.exit()
                    elif i > self.n_dipoles:
                        print('ERROR: change_signs changing sign to non-existing dipole number ' + str(i)) 
                        sys.exit()
                #
                signs_to_check = self.signs.copy()
                for index, dipole in enumerate(which_dipoles):
                    signs_to_check[dipole] = signs[0]
                self.check_and_assign_signs(signs_to_check)
        #
        else:
            print('ERROR: change_signs did not understand the input values')
            print('')
            print('HELP : which_dipoles = list of integers')
            print("       signs         = list of sign strings '+-' or '-+' or a single sign string")
            sys.exit()

            
    #
    # Check signs procedure
    #
    def check_and_assign_signs(self, signs = []):
        #
        """Procedure to make the sanity checks and assign the signs to a dipole object"""
        #
        if (type(signs) != list and type(signs) == str ):
            self.signs = [signs for i in range(0,self.n_dipoles)]
        #
        elif (type(signs) != list and type(signs) != str ):
            print('WARNING: check_and_assign_signs initializing a dipole object with an unvalid signs option. Usign default +-')
            self.signs = ['+-' for i in range(0,self.n_dipoles)]
        #
        elif (type(signs) == list and len(signs) == 0):
            self.signs = ['+-' for i in range(0,self.n_dipoles)]
        #
        elif (type(signs) == list and len(signs) == self.n_dipoles):
            #
            for sign in signs:
                #
                try:
                    value = sign.strip()
                except AttributeError:
                    print('ERROR: check_and_assign_signs initializing a dipole object with signs which are not understandable')
                    print('')
                    print('HELP : signs = list of sign strings (+- or -+) or a single sign strig')
                    sys.exit()
                #
                if (value != '+-' and value != '-+'):
                    print('ERROR: check_and_assign_signs initializing a dipole object with signs which are not understandable')
                    print('')
                    print('HELP : signs = list of sign strings (+- or -+) or a single sign strig')
                    sys.exit()
            #
            self.signs = [signs[i].strip() for i in range(0,self.n_dipoles)]
        #
        elif(type(signs) == list and len(signs) != self.n_dipoles and len(signs) != 0):
            print('ERROR: check_and_assign_signs initializing a dipole object with signs list shorter than the number of dipoles')
            sys.exit()
        #
        else:
            print('ERROR: check_and_assign_signs initializing a dipole object with not understandable signs')
            print('')
            print('HELP : signs = list of sign strings (+- or -+) or a single sign strig')
            sys.exit()
        #

    def position_the_dipoles_around(self, molecule, print_info = False):
        #
        """Procedure to generate the location of the fixed dipoles base on: """
        """ -The connectivity of the molecule
            -Its surface atoms
            -Other stuff"""
        #
        n_dipoles = 0
        positions = []
        versors   = []
        signs = '+-'
        #
        for surf_index in molecule.surface_atoms: 
            #
            # Sanity  check
            #
            if (molecule.number_connections[surf_index] == 0 and molecule.atomtypes[surf_index] != 'Cl'):
                print('ERROR: ' + atom + str(surf_index) + ' ' + str(molecule.atomtypes(surf_index)) + \
                      ' has zero connectivity.')
                sys.exit()
            #
            # If its a Cl atom, then generate the dipoles along the lone pairs directions
            #
            elif (molecule.number_connections[surf_index] == 0 and molecule.atomtypes[surf_index] == 'Cl'):
                #
                coords = molecule.coords[surf_index].copy()
                #
                A_coordinates = np.asarray([coords[0] + 1.0, coords[1] + 1.0, coords[2] + 1.0])
                B_coordinates = np.asarray([coords[0] - 1.0, coords[1] - 1.0, coords[2] + 1.0])
                C_coordinates = np.asarray([coords[0] - 1.0, coords[1] + 1.0, coords[2] - 1.0])
                D_coordinates = np.asarray([coords[0] + 1.0, coords[1] - 1.0, coords[2] - 1.0])
                #
                A_versor = np.asarray(molecule.coords[surf_index]     - A_coordinates)/\
                           np.linalg.norm(molecule.coords[surf_index] - A_coordinates)
                B_versor = np.asarray(molecule.coords[surf_index]     - B_coordinates)/\
                           np.linalg.norm(molecule.coords[surf_index] - B_coordinates)
                C_versor = np.asarray(molecule.coords[surf_index]     - C_coordinates)/\
                           np.linalg.norm(molecule.coords[surf_index] - C_coordinates)
                D_versor = np.asarray(molecule.coords[surf_index]     - D_coordinates)/\
                           np.linalg.norm(molecule.coords[surf_index] - D_coordinates)
                #
                positions.append(molecule.coords[surf_index] + A_versor)
                versors.append(A_versor)
                positions.append(molecule.coords[surf_index] + B_versor)
                versors.append(B_versor)
                positions.append(molecule.coords[surf_index] + C_versor)
                versors.append(C_versor)
                positions.append(molecule.coords[surf_index] + D_versor)
                versors.append(D_versor)
                n_dipoles += 4
                #
            #
            elif (molecule.number_connections[surf_index] == 1):
                #
                # Generate dipoles position:
                # 1) Along the line connceting it to the neighbours
                #
                connector_index = molecule.connected_to[surf_index][0]
                versor = np.asarray(molecule.coords[surf_index] - molecule.coords[connector_index])/\
                         np.linalg.norm(molecule.coords[surf_index] - molecule.coords[connector_index])
                #
                positions.append(molecule.coords[surf_index] + versor)
                versors.append(versor)
                n_dipoles += 1
                #
            elif (molecule.number_connections[surf_index] == 2):
                #
                # Generate dipoles position:
                # 1) Along the line connecting it to the neighbours
                # 2) Along the visector of the angle formed with the two neighbours
                # 3) Along the direction orthogonal to the plane identified by these atoms
                #
                connector_indices = molecule.connected_to[surf_index][:]
                #
                # Directions and in plane dipoles along connections
                #
                versor1 = np.asarray(molecule.coords[surf_index] - molecule.coords[connector_indices[0]])/\
                          np.linalg.norm(molecule.coords[surf_index] - molecule.coords[connector_indices[0]])
                #
                positions.append(molecule.coords[surf_index] + versor1)
                versors.append(versor1)
                n_dipoles += 1
                #
                versor2 = np.asarray(molecule.coords[surf_index] - molecule.coords[connector_indices[1]])/\
                          np.linalg.norm(molecule.coords[surf_index] - molecule.coords[connector_indices[1]])
                #
                positions.append(molecule.coords[surf_index] + versor2)
                versors.append(versor2)
                n_dipoles += 1
                #
                # Orthogonal dipoles
                #
                versor_ort = np.cross(versor1,versor2)
                versor_ort = versor_ort/np.linalg.norm(versor_ort)
                #
                positions.append(molecule.coords[surf_index] + versor_ort)
                versors.append(versor_ort)
                n_dipoles += 1
                #
                positions.append(molecule.coords[surf_index] - versor_ort)
                versors.append(-versor_ort)
                n_dipoles += 1
                #
                # In plane dipoles (bisectors)
                #
                versor3 = versor1 + versor2
                versor3 = versor3/np.linalg.norm(versor3)
                #
                positions.append(molecule.coords[surf_index] + versor3)
                versors.append(versor3)
                n_dipoles += 1
                #
                #positions.append(molecule.coords[surf_index] - versor3) #this is on the bisector but on the acute side
                # 
            elif (molecule.number_connections[surf_index] == 3):
                #
                # Generate dipoles position symmetrically over and below the plane identified by the atom and
                # the other atoms it is connected to. We try to use the backbone of the molecule to build the 
                # planes
                #
                connector_indices = molecule.connected_to[surf_index][:]
                connectors_not_on_surface = [i for i in connector_indices if i not in molecule.surface_atoms]
                #
                if (len(connectors_not_on_surface) != 2):
                    #
                    #Either 3 or less than 2 remaining neighbours: i take the first two to define the plane
                    #
                    connector_indices = connector_indices[0:2]
                else:
                    connector_indices = connectors_not_on_surface.copy()
                #
                # Now generate the direction orthogonal to the plane
                #
                versor1 = np.asarray(molecule.coords[surf_index] - molecule.coords[connector_indices[0]])/\
                          np.linalg.norm(molecule.coords[surf_index] - molecule.coords[connector_indices[0]])
                #
                versor2 = np.asarray(molecule.coords[surf_index] - molecule.coords[connector_indices[1]])/\
                          np.linalg.norm(molecule.coords[surf_index] - molecule.coords[connector_indices[1]])
                #
                versor3 = np.cross(versor1,versor2)
                versor3 = versor3/np.linalg.norm(versor3)
                #
                positions.append(molecule.coords[surf_index] + versor3)
                versors.append(versor3)
                n_dipoles += 1
                positions.append(molecule.coords[surf_index] - versor3)
                versors.append(-versor3)
                n_dipoles += 1
                #
            else:
                print("ERROR: more than three connecting atoms. Don't know what to do in this case.")
                sys.exit()
    
        self.n_dipoles = n_dipoles
        self.positions = positions.copy()
        self.directions = versors.copy()
        self.signs = [signs for i in range(0,n_dipoles)]
