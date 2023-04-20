from classes import molecule_class
from classes import dipoles_class
import constants
import numpy as np
import sys
import os
#
def create_EE_inp(et_seed_file = '', QMmolecule_file = '', EEdipoles_file = '', which_dipoles = [],\
                  et_output_file_name = '', target_directory =  './', computation_name = ''):
    #
    """Procedure to generate a .inp file (input for eT calculations)
       starting from a seed.inp file
       The seed.inp file, which is et_seed_file has everything set up except for:
       1) Name of the computation;
       2) Target QM molecule;
       3) Static dipole embedding the QM molecule; 
  
       The procedure is built in such a way that the user provides:
       1) et_seed_file.inp 
       2) QMmolecule_file.xyz (target QM molecule xyz file)
       3) EEdipoles_file.dip  (static dipolar embedding.dip file)
       4) et_output_file_name.inp (saved in target_directory)
       5) target_directory aforementioned (skip if you already have the path in the name of the file)
       6) which_dipoles (a list of dipoles you will include in the QM/EE calculation, default is all)
       7) computation_name (an ID inserted in the .inp file, useful for the user)


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
    et_create_EE_inp_sanity_checks(et_seed_file, QMmolecule_file, EEdipoles_file, \
                                   et_output_file_name, target_directory,     \
                                   computation_name, which_dipoles)
    #
    if target_directory[-1] != '/':
        target_directory +='/'
    #
    # Fetch the files
    #
    with open(et_seed_file, 'r') as seed:
        seed_lines = seed.readlines()
    #
    QMmolecule = molecule_class.molecule()
    QMmolecule.initialize_from_xyz(QMmolecule_file)
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
    # Initialize and work on the .inp file
    #
    with open(target_directory + et_output_file_name, 'w') as et_file:
        #
        # Fetch the basis used in the computation
        #
        basis  = 'aug-cc-pvdz'
        for line in seed_lines:
            if('basis:' in line):
                basis = line.split(':')[1]
        #
        # Write the .inp file with the correct information
        #
        for line in seed_lines:
            #
            if ('name:' in line):
                et_file.write('   name: ' + computation_name + '\n')
            #
            elif ('charge:' in line):
                et_file.write('   charge: ' + str(QMmolecule.charge) + '\n')
            #
            elif ('geometry' in line and 'end geometry' not in line):
                et_file.write('geometry \n')
                et_file.write('basis: '+ basis.strip() + '\n')
                for i, sym in enumerate(QMmolecule.atomtypes):
                    et_file.write(sym.rjust(2) + '  ' + \
                                  '{:5.5f}'.format(QMmolecule.coords[i][0]).rjust(10)+ '  ' + \
                                  '{:5.5f}'.format(QMmolecule.coords[i][1]).rjust(10)+ '  ' + \
                                  '{:5.5f}'.format(QMmolecule.coords[i][2]).rjust(10) + '\n')
                #
                et_file.write('--')
                #
                # EE embedding
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
                    et_file.write('\n' + 'H ' + ' [IMol=' + str(i+1).rjust(2) + ']  ' + \
                                   '{:5.5f}'.format(EEdipoles.positions[dip][0] + minus_vector[0]).rjust(10)+ '  ' + \
                                   '{:5.5f}'.format(EEdipoles.positions[dip][1] + minus_vector[1]).rjust(10)+ '  ' + \
                                   '{:5.5f}'.format(EEdipoles.positions[dip][2] + minus_vector[2]).rjust(10)+ '  ' + \
                                   '[q = ' + minus_vector_sign + '1.0] \n' )
                    #
                    # Charge in CM + distance
                    #
                    et_file.write( 'H ' + ' [IMol=' + str(i+1).rjust(2) + ']  ' + \
                                   '{:5.5f}'.format(EEdipoles.positions[dip][0] + plus_vector[0]).rjust(10)+ '  ' + \
                                   '{:5.5f}'.format(EEdipoles.positions[dip][1] + plus_vector[1]).rjust(10)+ '  ' + \
                                   '{:5.5f}'.format(EEdipoles.positions[dip][2] + plus_vector[2]).rjust(10)+ '  ' + \
                                   '[q = ' + plus_vector_sign + '1.0] \n' )
                    #
                    # Charge in CM, to be printed only as a debugging option (you need to uncomment it)
                    #
                    #et_file.write( 'H ' + ' [IMol=' + str(i+1).rjust(2) + ']  ' + \
                    #               '{:5.5f}'.format(EEdipoles.positions[dip][0]).rjust(10)+ '  ' + \
                    #               '{:5.5f}'.format(EEdipoles.positions[dip][1]).rjust(10)+ '  ' + \
                    #               '{:5.5f}'.format(EEdipoles.positions[dip][2]).rjust(10)+ '  ' + \
                    #               '[q = 0.0] \n' )
            #
            elif('xxx' in line):
                pass
            #
            elif('basis:' in line): #already done
                pass
            #
            else:
                et_file.write(line)
    #
#
#
#
def et_create_EE_inp_sanity_checks(et_seed_file, QMmolecule_file, EEdipoles_file, \
                                   et_output_file_name, target_directory,     \
                                   computation_name, which_dipoles):
    #
    """Routine for the sanity checks"""
    #
    # Check the types of the files
    #
    if (type(et_seed_file) != str):
        print('ERROR: create_EE_inp called with ununderstandable et_seed_file')
        sys.exit()
    #
    if (type(QMmolecule_file) != str):
        print('ERROR: create_EE_inp called with ununderstandable QMmolecule_file')
        sys.exit()
    #
    if (type(EEdipoles_file) != str):
        print('ERROR: create_EE_inp called with ununderstandable EEdipoles_file')
        sys.exit()
    #
    if (type(et_output_file_name) != str):
        print('ERROR: create_EE_inp called with ununderstandable output_file_name')
        sys.exit()
    #
    if (type(target_directory) != str):
        print('ERROR: create_EE_inp called with ununderstandable target_directory')
        sys.exit()
    #
    if (type(computation_name) != str):
        print('ERROR: create_EE_inp called with ununderstandable computation_name')
        sys.exit()
    #
    if (et_output_file_name.split('.')[-1] != 'inp'):
        print('ERROR: create_EE_inp called with output_file_name not ending with .inp')
        sys.exit()
    #
    # Check the which_dipoles type
    #
    elif(len(which_dipoles) != len(set(which_dipoles))):
        print('ERROR: create_EE_inp you have duplicates in which_dipoles' )
        sys.exit()
    #
    for i in which_dipoles:
        if type(i) != int:
            print('ERROR: create_EE_inp called with non-integer which_dipoles list')
            sys.exit()
    #
