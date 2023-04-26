from classes import molecule_class
from classes import dipoles_class
from classes import polarizable_embedding_class
import constants
import numpy as np
import sys
import os
#
def create_input(molecule_file = '', EEdipoles_file = '', which_dipoles = [], polarizable_embedding = [],\
                  nanofq_output_file_name = '', target_directory =  './', computation_name = ''):
    #
    """Procedure to generate a .mfq input file (input for nanoFQ calculations)
  
       The procedure is built in such a way that the user provides:
       1) molecule_file.xyz (target molecule xyz file treated at the FQ(FMu) level)
       2) EEdipoles_file.dip  (static dipolar embedding.dip file)
       3) nano_output_file_name.inp (saved in target_directory)
       4) target_directory aforementioned (skip if you already have the path in the name of the file)
       5) which_dipoles (a list of dipoles you will include in the FQ(FMu)/EE calculation, default is all)
       6) computation_name (an ID inserted in the .mfq file, useful for the user)
       7) polarizable_embedding (a polarizable_embedding class object with all the information for the computation)


       WARNING: all the files need to be provided wih the proper path to fetch them!
       WARNING: the .dip files contain information about the center of mass of the dipole
    """
    #
    distance = constants.dipoles_distance() #distance of the EE dipoles charges from their CM 
    #
    # Sanity checks
    #
    if type(which_dipoles) != list:
        which_dipoles = [which_dipoles]
    #
    nanofq_create_mfq_sanity_checks(molecule_file, EEdipoles_file,            \
                                    nanofq_output_file_name, target_directory,\
                                    computation_name, which_dipoles)
    #
    if target_directory[-1] != '/':
        target_directory +='/'
    #
    # Fetch the files
    #
    molecule = molecule_class.molecule()
    molecule.initialize_from_xyz(molecule_file)
    #
    EEdipoles = dipoles_class.dipoles()
    EEdipoles.initialize_from_dip(EEdipoles_file)
    #
    # Sanity checks on the initialized dipoles
    #
    if (len(which_dipoles) == 0): #displace them all
        which_dipoles = [i for i in range(0, EEdipoles.n_dipoles)]
    for i in which_dipoles:
        if i > EEdipoles.n_dipoles:
            print('ERROR: creat_EE_inp called with non-existing dipole number ' + str(i)) 
            sys.exit()
    #
    # Initialize the .mfq file and work on it
    #
    with open(target_directory + nanofq_output_file_name, 'w') as nano_file:
        #
        # Write the computation name
        #
        nano_file.write('!' + computation_name + '\n\n')
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
        nano_file.write('   static: ' + polarizable_embedding.force_field + '\n')
        nano_file.write('   kernel: gaussian\n') #TODO che kernel usare
        nano_file.write('end forcefield\n')
        #
        nano_file.write('\n')
        #
        # Write atomtypes section
        #
        nano_file.write('atom types\n')
        nano_file.write('   number: ' + str(len(polarizable_embedding.atomtypes)) + '\n')
        #
        for i,atomtype in enumerate(polarizable_embedding.atomtypes):
            string  = '   ' + atomtype + ': ['
            string += 'chi=' + str(polarizable_embedding.chi[i]) + ','
            string += 'eta=' + str(polarizable_embedding.eta[i])
            if (polarizable_embedding.force_field == 'fq'):
                string += ']\n'
            else: 
                string += ','
                if (polarizable_embedding.force_field == 'fq_pqeq'):
                    string += 'rq=' + str(polarizable_embedding.Rq[i]) + ']\n'
                elif(polarizable_embedding.force_field == 'fqfmu'):
                    string += 'alpha=' + str(polarizable_embedding.alpha[i]) + ']\n'
                elif(polarizable_embedding.force_field == 'fqfmu_pqeq'):
                    string += 'alpha=' + str(polarizable_embedding.alpha[i]) + ','
                    string += 'rq=' + str(polarizable_embedding.Rq[i]) + ','
                    string += 'rmu=' + str(polarizable_embedding.Rmu[i]) + ']\n'
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
        # Write the input geometry section
        #
        nano_file.write('input geometry \n')
        for i, sym in enumerate(molecule.atomtypes):
            nano_file.write(sym.rjust(2) + '  [IMol=1] ' + \
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
        for i,dip in enumerate(which_dipoles):
            # 
            plus_vector  =  distance*EEdipoles.directions[dip,:]
            minus_vector = -distance*EEdipoles.directions[dip,:]
            plus_vector_sign  = EEdipoles.signs[dip][1]
            minus_vector_sign = EEdipoles.signs[dip][0]
            #
            # Charge in CM - distance
            #
            nano_file.write('XX' + \
                            '{:5.5f}'.format(EEdipoles.positions[dip][0] + minus_vector[0]).rjust(10)+ '  ' + \
                            '{:5.5f}'.format(EEdipoles.positions[dip][1] + minus_vector[1]).rjust(10)+ '  ' + \
                            '{:5.5f}'.format(EEdipoles.positions[dip][2] + minus_vector[2]).rjust(10)+ '  ' + \
                            minus_vector_sign + '1.0\n' )
            #
            # Charge in CM + distance
            #
            nano_file.write('XX' + \
                            '{:5.5f}'.format(EEdipoles.positions[dip][0] + plus_vector[0]).rjust(10)+ '  ' + \
                            '{:5.5f}'.format(EEdipoles.positions[dip][1] + plus_vector[1]).rjust(10)+ '  ' + \
                            '{:5.5f}'.format(EEdipoles.positions[dip][2] + plus_vector[2]).rjust(10)+ '  ' + \
                            plus_vector_sign + '1.0\n' )
            #
            # Charge in CM, to be printed only as a debugging option (you need to uncomment it)
            #
            #nano_file.write( 'H ' + ' [IMol=' + str(i+1).rjust(2) + ']  ' + \
            #               '{:5.5f}'.format(EEdipoles.positions[dip][0]).rjust(10)+ '  ' + \
            #               '{:5.5f}'.format(EEdipoles.positions[dip][1]).rjust(10)+ '  ' + \
            #               '{:5.5f}'.format(EEdipoles.positions[dip][2]).rjust(10)+ '  ' + \
            #               '[q = 0.0] \n' )
            #
        #
        nano_file.write('end electrostatic embedding geometry\n')
        #
#
#
#
def nanofq_create_mfq_sanity_checks(molecule_file, EEdipoles_file, \
                                   nanofq_output_file_name, target_directory,     \
                                   computation_name, which_dipoles):
    #
    """Routine for the sanity checks"""
    #
    # Check the types of the files
    #
    if (type(molecule_file) != str):
        print('ERROR: create_mfq called with ununderstandable QMmolecule_file')
        sys.exit()
    #
    if (type(EEdipoles_file) != str):
        print('ERROR: create_mfq called with ununderstandable EEdipoles_file')
        sys.exit()
    #
    if (type(nanofq_output_file_name) != str):
        print('ERROR: create_mfq called with ununderstandable output_file_name')
        sys.exit()
    #
    if (type(target_directory) != str):
        print('ERROR: create_mfq called with ununderstandable target_directory')
        sys.exit()
    #
    if (type(computation_name) != str):
        print('ERROR: create_mfq called with ununderstandable computation_name')
        sys.exit()
    #
    if (nanofq_output_file_name.split('.')[-1] != 'mfq'):
        print('ERROR: create_mfq called with output_file_name not ending with .mfq')
        sys.exit()
    #
    # Check the which_dipoles type
    #
    elif(len(which_dipoles) != len(set(which_dipoles))):
        print('ERROR: create_mfq you have duplicates in which_dipoles' )
        sys.exit()
    #
    for i in which_dipoles:
        if type(i) != int:
            print('ERROR: create_mfq called with non-integer which_dipoles list')
            sys.exit()
    #
#
# Get_energy procedure
#
def get_energy(log_file = ''):
    #
    """Procedure to read a nanofq output .log file and get the interaction energy FQ(FMu)/EE"""
    #
    # Fetch and open file
    #
    if( not log_file.endswith('.log')):
        print('ERROR: get_energy without a .log file provided')
        sys.exit()
    #
    with open(log_file, 'r') as file_:
        lines = file_.readlines()
    #
    for line in lines:
        if ('Electrostatic Embedding Interaction =') in line:
            energy = float(line.split('=')[1]) 
            break
    #
    return energy
        
