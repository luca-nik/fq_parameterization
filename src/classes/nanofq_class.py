import numpy as np
import constants
from classes import molecule_class
from classes import dipoles_class
from classes import cluster_class
from classes import polarizable_embedding_class
import sys
import subprocess
#
class nanofq:
    #
    """nanofq class object"""
    #
    def __init__(self, nanofq_path = '',\
                 polarizable_model = polarizable_embedding_class.polarizable_embedding(),\
                 molecule = molecule_class.molecule(), dipoles = dipoles_class.dipoles()):
        #
        """Initialization of the nanofq_class object
            nanofq_path       = str (the path where you can find the nanoFQ code)
            polarizable_model = polarizable_embedding_class object
            molecule          = molecule_class object 
            dipoles           = dipoles_class object 
        """
        #
        # Sanity checks and initialization
        #
        self.check_inputs(nanofq_path,polarizable_model, molecule, dipoles)
        #
        self.input   = '' #.mfq
        self.output  = '' #.log
        self.comment = ''
        self.name    = ''
        #
        self.which_dipoles = [] #which dipoles to include in the calculation
        #
    #
    def print_info(self):
        #
        print('***********************************************************')
        print('nanoFQ object: ')
        print('path          : ' + self.nanofq_path)
        print('')
        if isinstance(self.molecule, molecule_class.molecule):
           print('System        : molecule')
           print('-----------------------')
           self.molecule.print_info()
        #
        elif isinstance(self.molecule, cluster_class.cluster):
           print('System        : cluster')
           print('-----------------------')
           for mol in self.molecule.molecules:
               mol.print_info()
               print('-----------------------')
        print('')
        if (self.dipoles.n_dipoles != 0):
            if (type(self.which_dipoles) == list and len(self.which_dipoles) != 0):
                print('EE dipoles      : ' + str(len(self.which_dipoles)))
                print('-----------------------')
                print('Dipole | ' + 'Pos X'.rjust(10) + ' ' +'Pos Y'.rjust(10) + ' ' + 'Pos Z'.rjust(10) + \
                               ' | ' + 'Dir X'.rjust(10) + ' ' +'Dir Y'.rjust(10) + ' ' + 'Dir Z'.rjust(10)      + \
                               ' | ' + 'Sign'.rjust(6))
                for i in self.which_dipoles:
                    print(str(i).ljust(6) + ' | ' + \
                             '{:5.5f}'.format(self.dipoles.positions[i][0]).rjust(10)+ ' '    + \
                             '{:5.5f}'.format(self.dipoles.positions[i][1]).rjust(10)+ ' '    + \
                             '{:5.5f}'.format(self.dipoles.positions[i][2]).rjust(10) + ' | ' + \
                             '{:5.5f}'.format(self.dipoles.directions[i][0]).rjust(10)+ ' '   + \
                             '{:5.5f}'.format(self.dipoles.directions[i][1]).rjust(10)+ ' '   + \
                             '{:5.5f}'.format(self.dipoles.directions[i][2]).rjust(10)+ ' | ' + \
                             self.dipoles.signs[i].rjust(6)
                             )
            elif (type(self.which_dipoles) == int):
                print('EE dipoles      : 1')
                print('-----------------------')
                print('Dipole | ' + 'Pos X'.rjust(10) + ' ' +'Pos Y'.rjust(10) + ' ' + 'Pos Z'.rjust(10) + \
                               ' | ' + 'Dir X'.rjust(10) + ' ' +'Dir Y'.rjust(10) + ' ' + 'Dir Z'.rjust(10)      + \
                               ' | ' + 'Sign'.rjust(6))
                #
                i = self.which_dipoles
                print(str(i).ljust(6) + ' | ' + \
                         '{:5.5f}'.format(self.dipoles.positions[i][0]).rjust(10)+ ' '    + \
                         '{:5.5f}'.format(self.dipoles.positions[i][1]).rjust(10)+ ' '    + \
                         '{:5.5f}'.format(self.dipoles.positions[i][2]).rjust(10) + ' | ' + \
                         '{:5.5f}'.format(self.dipoles.directions[i][0]).rjust(10)+ ' '   + \
                         '{:5.5f}'.format(self.dipoles.directions[i][1]).rjust(10)+ ' '   + \
                         '{:5.5f}'.format(self.dipoles.directions[i][2]).rjust(10)+ ' | ' + \
                         self.dipoles.signs[i].rjust(6)
                         )
            #
            #
            #
            else:
                print('EE dipoles      : ' + str(self.dipoles.n_dipoles))
                print('-----------------------')
                print('Dipole | ' + 'Pos X'.rjust(10) + ' ' +'Pos Y'.rjust(10) + ' ' + 'Pos Z'.rjust(10) + \
                               ' | ' + 'Dir X'.rjust(10) + ' ' +'Dir Y'.rjust(10) + ' ' + 'Dir Z'.rjust(10)      + \
                               ' | ' + 'Sign'.rjust(6))
                for i in range(self.dipoles.n_dipoles):
                    print(str(i).ljust(6) + ' | ' + \
                             '{:5.5f}'.format(self.dipoles.positions[i][0]).rjust(10)+ ' '    + \
                             '{:5.5f}'.format(self.dipoles.positions[i][1]).rjust(10)+ ' '    + \
                             '{:5.5f}'.format(self.dipoles.positions[i][2]).rjust(10) + ' | ' + \
                             '{:5.5f}'.format(self.dipoles.directions[i][0]).rjust(10)+ ' '   + \
                             '{:5.5f}'.format(self.dipoles.directions[i][1]).rjust(10)+ ' '   + \
                             '{:5.5f}'.format(self.dipoles.directions[i][2]).rjust(10)+ ' | ' + \
                             self.dipoles.signs[i].rjust(6)
                             )
        
        print('')
        print('PE           :')
        self.polarizable_model.print_info()
        print('***********************************************************')
        #
        #
        #
        pass
#
    def create_ee_input(self,input_ = '', computation_comment = '', which_dipoles = []):
        #
        """Procedure to generate a .mfq input file (input for nanoFQ calculations)
           1) input_ the complete path + name of the .mfq file you want to create
           2) computation_name (an ID inserted in the .mfq file, useful for the user)
           3) which_dipoles (a list of dipoles you will include in the FQ(FMu)/EE calculation, default is all)
        """
        #
        distance = constants.dipoles_distance() #distance of the EE dipoles charges from their CM 
        #
        # Sanity checks
        #
        if type(which_dipoles) != list:
            which_dipoles = [which_dipoles]
        #
        self.mfq_sanity_checks(input_,computation_comment, which_dipoles)
        #
        # Sanity checks on the initialized dipoles
        #
        if (len(which_dipoles) == 0): #displace them all
            which_dipoles = [i for i in range(0, self.dipoles.n_dipoles)]
        for i in which_dipoles:
            if i > self.dipoles.n_dipoles:
                print('ERROR: create_ee_input called with non-existing dipole number ' + str(i)) 
                sys.exit()
        #
        # Initialize the .mfq file and work on it
        #
        with open(self.input, 'w') as nano_file:
            #
            # Write the computation comment
            #
            nano_file.write('!' + self.comment + '\n\n')
            #
            # Write what section
            #
            nano_file.write('what\n')
            nano_file.write('   energy\n')
            nano_file.write('end what\n')
            #
            nano_file.write('\n')
            #
            # Write force field section
            #
            nano_file.write('forcefield\n')
            nano_file.write('   static: ' + self.polarizable_model.force_field + '\n')
            if (self.polarizable_model.force_field == 'fq' or self.polarizable_model.force_field == 'fq_pqeq'):
                nano_file.write('   kernel: ohno\n') 
            else:
                nano_file.write('   kernel: gaussian\n') 

            nano_file.write('end forcefield\n')
            #
            nano_file.write('\n')
            #
            # Write atomtypes section
            #
            nano_file.write('atom types\n')
            nano_file.write('   number: ' + str(len(self.polarizable_model.atomtypes)) + '\n')
            #
            for i,atomtype in enumerate(self.polarizable_model.atomtypes):
                string  = '   ' + atomtype + ': ['
                string += 'chi=' + str(self.polarizable_model.chi[i]) + ','
                string += 'eta=' + str(self.polarizable_model.eta[i])
                if (self.polarizable_model.force_field == 'fq'):
                    string += ']\n'
                else: 
                    string += ','
                    if (self.polarizable_model.force_field == 'fq_pqeq'):
                        string += 'rq=' + str(self.polarizable_model.Rq[i]) + ']\n'
                    elif(self.polarizable_model.force_field == 'fqfmu'):
                        string += 'alpha=' + str(self.polarizable_model.alpha[i]) + ']\n'
                    elif(self.polarizable_model.force_field == 'fqfmu_pqeq'):
                        string += 'alpha=' + str(self.polarizable_model.alpha[i]) + ','
                        string += 'rq=' + str(self.polarizable_model.Rq[i]) + ','
                        string += 'rmu=' + str(self.polarizable_model.Rmu[i]) + ']\n'
                #
                nano_file.write(string)
            #
            nano_file.write('end atom types\n')
            #
            nano_file.write('\n')
            #
            # Write output section
            #
            nano_file.write('output\n')
            nano_file.write('   verbose: 1\n')
            nano_file.write('end output\n')
            #
            nano_file.write('\n')
            #
            # Write the input geometry section paying attention wether is a list of molecules or just one
            #
            nano_file.write('input geometry \n')
            #
            if (not isinstance(self.molecule,cluster_class.cluster)):
                for i, sym in enumerate(self.molecule.atomtypes):
                    nano_file.write(sym.rjust(2) + '  [IMol=1] ' + \
                                  '{:5.5f}'.format(self.molecule.coords[i][0]).rjust(10)+ '  ' + \
                                  '{:5.5f}'.format(self.molecule.coords[i][1]).rjust(10)+ '  ' + \
                                  '{:5.5f}'.format(self.molecule.coords[i][2]).rjust(10) + '\n')
            else:
                for mol_indx, molecule in enumerate(self.molecule.molecules):
                    for i, sym in enumerate(molecule.atomtypes):
                        if (molecule.charge == 0):
                            nano_file.write(sym.rjust(2) + '  [IMol=' + str(mol_indx + 1).ljust(2) + '] '  + \
                                          '{:5.5f}'.format(molecule.coords[i][0]).rjust(10)+ '  ' + \
                                          '{:5.5f}'.format(molecule.coords[i][1]).rjust(10)+ '  ' + \
                                          '{:5.5f}'.format(molecule.coords[i][2]).rjust(10) + '\n')
                        else:
                            nano_file.write(sym.rjust(2) + '  [IMol=' + str(mol_indx + 1).ljust(2) + \
                                            '|charge=' + str(molecule.charge).ljust(4) + '] '  + \
                                          '{:5.5f}'.format(molecule.coords[i][0]).rjust(10)+ '  ' + \
                                          '{:5.5f}'.format(molecule.coords[i][1]).rjust(10)+ '  ' + \
                                          '{:5.5f}'.format(molecule.coords[i][2]).rjust(10) + '\n')
            #
            nano_file.write('end input geometry\n')
            #
            nano_file.write('\n')
            #
            # Write the electostatic embedding geometry
            #
            nano_file.write('electrostatic embedding geometry\n')
            #
            for i,dip in enumerate(self.which_dipoles):
                # 
                plus_vector  =  distance*self.dipoles.directions[dip,:]
                minus_vector = -distance*self.dipoles.directions[dip,:]
                plus_vector_sign  = self.dipoles.signs[dip][1]
                minus_vector_sign = self.dipoles.signs[dip][0]
                #
                # Charge in CM - distance
                #
                nano_file.write('XX' + \
                                '{:5.5f}'.format(self.dipoles.positions[dip][0] + minus_vector[0]).rjust(10)+ '  ' + \
                                '{:5.5f}'.format(self.dipoles.positions[dip][1] + minus_vector[1]).rjust(10)+ '  ' + \
                                '{:5.5f}'.format(self.dipoles.positions[dip][2] + minus_vector[2]).rjust(10)+ '  ' + \
                                minus_vector_sign + '1.0\n' )
                #
                # Charge in CM + distance
                #
                nano_file.write('XX' + \
                                '{:5.5f}'.format(self.dipoles.positions[dip][0] + plus_vector[0]).rjust(10)+ '  ' + \
                                '{:5.5f}'.format(self.dipoles.positions[dip][1] + plus_vector[1]).rjust(10)+ '  ' + \
                                '{:5.5f}'.format(self.dipoles.positions[dip][2] + plus_vector[2]).rjust(10)+ '  ' + \
                                plus_vector_sign + '1.0\n' )
                #
                # Charge in CM, to be printed only as a debugging option (you need to uncomment it)
                #
                #nano_file.write( 'H ' + ' [IMol=' + str(i+1).rjust(2) + ']  ' + \
                #               '{:5.5f}'.format(self.dipoles.positions[dip][0]).rjust(10)+ '  ' + \
                #               '{:5.5f}'.format(self.dipoles.positions[dip][1]).rjust(10)+ '  ' + \
                #               '{:5.5f}'.format(self.dipoles.positions[dip][2]).rjust(10)+ '  ' + \
                #               '[q = 0.0] \n' )
                #
            #
            nano_file.write('end electrostatic embedding geometry\n')
            #
#
    def create_polar_input(self,input_ = '', computation_comment = ''):
        #
        """Procedure to generate a .mfq input file for the computation of the polar(input for nanoFQ calculations)
           1) input_ the complete path + name of the .mfq file you want to create
           2) computation_name (an ID inserted in the .mfq file, useful for the user)
        """
        #
        # Sanity checks
        #
        self.mfq_sanity_checks(input_,computation_comment)
        #
        # Initialize the .mfq file and work on it
        #
        with open(self.input, 'w') as nano_file:
            #
            # Write the computation comment
            #
            nano_file.write('!' + self.comment + '\n\n')
            #
            # Write what section
            #
            nano_file.write('what\n')
            nano_file.write('   static response\n')
            nano_file.write('end what\n')
            #
            nano_file.write('\n')
            #
            # Write force field section
            #
            nano_file.write('forcefield\n')
            nano_file.write('   static: ' + self.polarizable_model.force_field + '\n')
            if (self.polarizable_model.force_field == 'fq' or self.polarizable_model.force_field == 'fq_pqeq'):
                nano_file.write('   kernel: ohno\n') 
            else:
                nano_file.write('   kernel: gaussian\n') 
            nano_file.write('end forcefield\n')
            #
            nano_file.write('\n')
            #
            # Write atomtypes section
            #
            nano_file.write('atom types\n')
            nano_file.write('   number: ' + str(len(self.polarizable_model.atomtypes)) + '\n')
            #
            for i,atomtype in enumerate(self.polarizable_model.atomtypes):
                string  = '   ' + atomtype + ': ['
                string += 'chi=' + str(self.polarizable_model.chi[i]) + ','
                string += 'eta=' + str(self.polarizable_model.eta[i])
                if (self.polarizable_model.force_field == 'fq'):
                    string += ']\n'
                else: 
                    string += ','
                    if (self.polarizable_model.force_field == 'fq_pqeq'):
                        string += 'rq=' + str(self.polarizable_model.Rq[i]) + ']\n'
                    elif(self.polarizable_model.force_field == 'fqfmu'):
                        string += 'alpha=' + str(self.polarizable_model.alpha[i]) + ']\n'
                    elif(self.polarizable_model.force_field == 'fqfmu_pqeq'):
                        string += 'alpha=' + str(self.polarizable_model.alpha[i]) + ','
                        string += 'rq=' + str(self.polarizable_model.Rq[i]) + ','
                        string += 'rmu=' + str(self.polarizable_model.Rmu[i]) + ']\n'
                #
                nano_file.write(string)
            #
            nano_file.write('end atom types\n')
            #
            nano_file.write('\n')
            #
            # Write the field section
            #
            nano_file.write('field\n')
            nano_file.write('   efield: static\n')
            nano_file.write('   rhs: potential\n')
            nano_file.write('   field intensity: 0.1\n')
            nano_file.write('end field\n')
            #
            nano_file.write('\n')
            #
            # Write output section
            #
            nano_file.write('output\n')
            nano_file.write('   verbose: 1\n')
            nano_file.write('end output\n')
            #
            nano_file.write('\n')
            #
            # Write the input geometry section paying attention wether is a list of molecules or just one
            #
            nano_file.write('input geometry \n')
            #
            if (not isinstance(self.molecule, cluster_class.cluster)):
                for i, sym in enumerate(self.molecule.atomtypes):
                    nano_file.write(sym.rjust(2) + '  [IMol=1] ' + \
                                  '{:5.5f}'.format(self.molecule.coords[i][0]).rjust(10)+ '  ' + \
                                  '{:5.5f}'.format(self.molecule.coords[i][1]).rjust(10)+ '  ' + \
                                  '{:5.5f}'.format(self.molecule.coords[i][2]).rjust(10) + '\n')
            else:
                for mol_indx, molecule in enumerate(self.molecule.molecules):
                    for i, sym in enumerate(molecule.atomtypes):
                        nano_file.write(sym.rjust(2) + '  [IMol=' + str(mol_indx + 1).ljust(2) + '] '  + \
                                      '{:5.5f}'.format(molecule.coords[i][0]).rjust(10)+ '  ' + \
                                      '{:5.5f}'.format(molecule.coords[i][1]).rjust(10)+ '  ' + \
                                      '{:5.5f}'.format(molecule.coords[i][2]).rjust(10) + '\n')
            #
            nano_file.write('end input geometry\n')
            #
    #
    # Procedure for the input sanity checks
    #
    def mfq_sanity_checks(self,input_,computation_comment, which_dipoles=[]):
        #
        """Routine for the sanity checks"""
        #
        if (type(input_) != str):
            print('ERROR: create_mfq called with ununderstandable input_ name')
            sys.exit()
        #
        if (type(computation_comment) != str):
            print('ERROR: create_mfq called with ununderstandable computation_comment')
            sys.exit()
        #
        if (input_.split('.')[-1] == input_):
            input_ = input_ + '.mfq'
        elif (input_.split('.')[-1] != 'mfq' and input_.split('.')[-1] != input_):
            print('ERROR: create_mfq called with weird input. Avoid "input_ = file.stuff" and use "input_ = file.mfq" or "input_ = file"')
            sys.exit()
        #
        # Check the which_dipoles type
        #
        elif(len(which_dipoles) != len(set(which_dipoles))):
            print('ERROR: create_mfq you have duplicates in which_dipoles' )
            sys.exit()
        #
        if which_dipoles != []:
            for i in which_dipoles:
                if type(i) != int:
                    print('ERROR: create_mfq called with non-integer which_dipoles list')
                    sys.exit()
        #
        self.input = input_
        self.which_dipoles = which_dipoles
        self.comment = computation_comment
    #
    # Procedure to run the nanofq_code
    #
    def run(self):
        subprocess.run([self.nanofq_path, self.input])
        if self.name == '':
            self.guess_name_from_input()
        self.output = self.name + '.log'
    #
    # Guess the name of the file
    #
    def guess_name_from_input(self):
        self.name = self.input.split('.mfq')[0]
    #
    # Guess the name of the file
    #
    def guess_name_from_dip(self):
        nomignolo = self.dipoles.name.split('.dip')[0]
        self.name = nomignolo.split('/')[-1]
        return self.name
    #
    # Get_energy procedure
    #
    def get_energy(self, which = 'electrostatic'):
        #
        """Procedure to read a nanofq output .log file and get the interaction energy FQ(FMu)/EE"""
        #
        # Fetch and open file
        #
        if( not self.output.endswith('.log')):
            print('ERROR: get_energy without a .log file provided')
            sys.exit()
        #
        with open(self.output, 'r') as file_:
            lines = file_.readlines()
        #
        if (which.lower() == 'electrostatic'):
            target_string = 'Electrostatic Embedding Interaction ='
        elif (which.lower() == 'total'):
            target_string = 'Energy ='
        else:
            print('ERROR: get_energy without proper energy requested')
            print('       Allowed keywords: electrostatic, total ')
            sys.exit()
        #
        for line in lines:
            if (target_string) in line:
                energy_line = line.split('=')[1]
                energy = float(energy_line.split('a.u.')[0])
                return energy 
        print('ERROR: no energy found')
        sys.exit()
    #
    # Get_polar procedure
    #
    def get_polar(self, which = 'isotropic'):
        #
        """Procedure to read a nanofq output .log file and get the polarizability"""
        #
        # Fetch and open file
        #
        if( not self.output.endswith('.log')):
            print('ERROR: get_polar without a .log file provided')
            sys.exit()
        #
        with open(self.output, 'r') as file_:
            lines = file_.readlines()
        #
        if (which.lower() == 'isotropic'):
            target_string = 'Polar Iso  ='
        elif (which.lower() == 'tensor'):
            target_string = 'Polarizability Tensor'
        else:
            print('ERROR: get_polar without proper polarizability requested.')
            print('       Allowed keywords: isotropic, tensor ')
            sys.exit()
        #
        for indx,line in enumerate(lines):
            if (target_string == 'Polar Iso  =' and target_string in line):
                polar_line = line.split('=')[1]
                polar = float(polar_line.split('a.u.')[0])
                return polar 
            elif(target_string == 'Polarizability Tensor' and target_string in line):
                polar = []
                x_values = [float(i) for i in lines[indx+3].split()[1:]]
                y_values = [float(i) for i in lines[indx+4].split()[1:]]
                z_values = [float(i) for i in lines[indx+5].split()[1:]]
                polar.append(x_values)
                polar.append(y_values)
                polar.append(z_values)
                return np.asarray(polar)
                
        print('ERROR: no polar found')
        sys.exit()
    #
    #
    #
    def check_inputs(self,nanofq_path,polarizable_model,molecule,dipoles):
        #
        """Routine for the sanity checks"""
        #
        # Check the types of the input objects
        #
        if (type(nanofq_path) != str):
            print('ERROR: nanofq_class initialization without proper nanofq_path')
            sys.exit()
        else:
            self.nanofq_path = nanofq_path
        #
        # Polarizable model
        #
        if (not isinstance(polarizable_model, polarizable_embedding_class.polarizable_embedding)):
            print('ERROR: nanofq_class initialization without proper polarizable model specified') 
            sys.exit()
        else:
            self.polarizable_model = polarizable_model
        #
        if (isinstance(molecule, cluster_class.cluster)):
            for mol in molecule.molecules:
                if (not isinstance(mol, molecule_class.molecule)):
                    print('ERROR: nanofq_class initialization without proper molecule object specified') 
                    sys.exit()
            else:
                self.molecule = molecule
        #
        elif (not isinstance(molecule, molecule_class.molecule) and not isinstance(molecule,cluster_class.cluster)):
            print('ERROR: nanofq_class initialization without proper molecule object specified') 
            sys.exit()
        else:
            self.molecule = molecule
        #
        if (not isinstance(dipoles, dipoles_class.dipoles)):
            print('ERROR: nanofq_class initialization without proper dipoles object specified') 
            sys.exit()
        else:
            self.dipoles = dipoles
