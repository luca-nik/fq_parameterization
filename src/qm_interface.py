from classes import molecule_class
from classes import dipoles_class
from classes import cluster_class
import numpy as np
import sys
import os
#
def create_et_EE_inp(et_seed_file = '', QMmolecule = '', EEdipoles = '', which_dipoles = [],\
                     et_output_file_name = '', target_directory =  './', computation_name = '',\
                     dipoles_distance = 0.5, single_charge = [False, 0.0]):
    #
    """Procedure to generate a .inp file (input for eT calculations)
       starting from a seed.inp file
       The seed.inp file, which is et_seed_file has everything set up except for:
       1) Name of the computation;
       2) Target QM molecule;
       3) Static dipole embedding of the QM molecule; 
  
       The procedure is built in such a way that the user provides:
         1) et_seed_file.inp        : string (This is a file where we read all the QM computation details)
         2) QMmolecule              : molecule class (Target QM molecule)
         3) EEdipoles               : dipoles_class  (Electrostatic embeddign dipoles)
         4) et_output_file_name.inp : string (output .inp file name, saved in target_directory)
         5) target_directory        : string (where you save the .inp files,
                                              skip if you already have the path in the name of the file)
         6) which_dipoles           : list of int (list of the indexes of the dipoles to put in the 
                                                   QM/EE calculation, default is all)
         7) computation_name        : string  (an ID inserted in the .inp file, useful for the user)
         8) dipoles_distance        : float (distance of each charge of the dipole from the dipole
                                             center of mass)

       WARNING: all the files need to be provided wih the proper path to fetch them!
    """
    #
    distance = dipoles_distance
    #
    # Sanity checks
    #
    if type(which_dipoles) != list:
        which_dipoles = [which_dipoles]
    #
    et_create_et_EE_inp_sanity_checks(et_seed_file, QMmolecule, EEdipoles, \
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
                et_file.write('   name: ' + computation_name + ' \n')
            #
            elif ('charge:' in line):
                et_file.write('   charge: ' + str(int(QMmolecule.charge)) + '\n')
            #
            elif ('geometry' in line and 'end geometry' not in line):
                et_file.write('geometry \n')
                et_file.write('basis: '+ basis.strip() + '\n')
                #
                # QM molecule
                #
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
                if single_charge[0]:
                    for i,dip in enumerate(which_dipoles):
                        #
                        # Charge in CM
                        #
                        et_file.write( 'H ' + ' [IMol=' + str(i+1).rjust(2) + ']  ' + \
                                       '{:5.5f}'.format(EEdipoles.positions[dip][0]).rjust(10)+ '  ' + \
                                       '{:5.5f}'.format(EEdipoles.positions[dip][1]).rjust(10)+ '  ' + \
                                       '{:5.5f}'.format(EEdipoles.positions[dip][2]).rjust(10)+ '  ' + \
                                       '[q = ' + str(single_charge[1]) + '] \n' )
                else:
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
def et_create_et_EE_inp_sanity_checks(et_seed_file, QMmolecule, EEdipoles, \
                                   et_output_file_name, target_directory,     \
                                   computation_name, which_dipoles):
    #
    """Routine for the sanity checks"""
    #
    # Check the types of the files
    #
    if (type(et_seed_file) != str):
        print('ERROR: create_et_EE_inp called with ununderstandable et_seed_file')
        sys.exit()
    #
    if (not isinstance(QMmolecule,molecule_class.molecule)):
        print('ERROR: create_et_EE_inp called without proper QM molecule object')
        sys.exit()
    #
    if (not isinstance(EEdipoles,dipoles_class.dipoles)):
        print('ERROR: create_et_EE_inp called without proper EE dipoles object')
        sys.exit()
    #
    if (type(et_output_file_name) != str):
        print('ERROR: create_et_EE_inp called with ununderstandable output_file_name')
        sys.exit()
    #
    if (type(target_directory) != str):
        print('ERROR: create_et_EE_inp called with ununderstandable target_directory')
        sys.exit()
    #
    if (type(computation_name) != str):
        print('ERROR: create_et_EE_inp called with ununderstandable computation_name')
        sys.exit()
    #
    if (et_output_file_name.split('.')[-1] != 'inp'):
        print('ERROR: create_et_EE_inp called with output_file_name not ending with .inp')
        sys.exit()
    #
    # Check the which_dipoles type
    #
    elif(len(which_dipoles) != len(set(which_dipoles))):
        print('ERROR: create_et_EE_inp you have duplicates in which_dipoles' )
        sys.exit()
    #
    for i in which_dipoles:
        if type(i) != int:
            print('ERROR: create_et_EE_inp called with non-integer which_dipoles list')
            sys.exit()
    #
#
#
#
def et_read_EE_energy(et_out_file = ''):
    #
    """Routine to read the electrostatic energy of et """
    #
    #
    #
    try:
        if (et_out_file.split('.')[-1] != 'out'):
            print('ERROR: et_read_EE_energy witout proper .out file provided')
            sys.exit()
    except:
        print('ERROR: et_read_EE_energy witout proper .out file provided')
        sys.exit()
    #
    with open(et_out_file, 'r') as out:
        lines = out.readlines()
    #
    for line in lines:
        if ('QM/MM Electrostatic Energy') in line:
            linea = line.split('QM/MM Electrostatic Energy:')[-1]
            energy = float(linea.split()[0])
            break
    # 
    return energy
#
#
#
def adf_read_polar(adf_out_file = '', which = 'isotropic'):
    #
    """Routine to read the polarizability from adf file """
    #
    #
    #
    try:
        if (adf_out_file.split('.')[-1] != 'log'):
            print('ERROR: adf_read_polar witout proper .log file provided')
            sys.exit()
    except:
        print('ERROR: adf_read_polar witout proper .log file provided')
        sys.exit()
    #
    with open(adf_out_file, 'r') as out:
        lines = out.readlines()
    #
    targets_list = ['Isotropic        POLARIZABILITY   =', 'Polarizability tensor:']
    if which == 'isotropic':
        target = targets_list[0]
    elif which == 'tensor':
        target = targets_list[1] 
    else:
        print('ERROR: adf_read_polar without proper polar requested. Possibilities isotropic,tensor')
    #
    for indx,line in enumerate(lines):
        if target in line:
            if which == 'isotropic':
                polar = float(line.split(target)[1].split()[0])
                break
            else:
                polar = []
                x = [float(i) for i in lines[indx+2].split()]
                y = [float(i) for i in lines[indx+3].split()]
                z = [float(i) for i in lines[indx+4].split()]
                polar.append(x)
                polar.append(y)
                polar.append(z)
                polar = np.asarray(polar)
                break
    # 
    return polar
#
#
#
def create_adf_polar_inp(adf_seed_file = '', cluster = '', adf_output_file_name = '', target_directory =  './'):
    #
    # Sanity checks
    #
    create_adf_polar_inp_sanity_checks(adf_seed_file,cluster,adf_output_file_name,target_directory)
    #
    if target_directory[-1] != '/':
        target_directory +='/'
    #
    # Fetch the files
    #
    with open(adf_seed_file, 'r') as seed:
        seed_lines = seed.readlines()
    #
    # Initialize and work on the .inp file
    #
    with open(target_directory + adf_output_file_name, 'w') as adf_file:
        #
        # Write down information from seed file
        #
        for line in seed_lines:
            if('  Atoms' in line):
                adf_file.write(line)
                for mol_indx, molecule in enumerate(cluster.molecules):
                    for i, sym in enumerate(molecule.atomtypes):
                        #
                        # Take the right symbol
                        #
                        if sym == 'Cl':
                            right_sym = sym
                        else:
                            right_sym = sym[0]
                        #
                        adf_file.write('    ' + right_sym.ljust(2)  + \
                                       '{:5.5f}'.format(molecule.coords[i][0]).rjust(10)+ '  ' + \
                                       '{:5.5f}'.format(molecule.coords[i][1]).rjust(10)+ '  ' + \
                                       '{:5.5f}'.format(molecule.coords[i][2]).rjust(10) + '\n')
            else:
                adf_file.write(line)
#
#
#
def create_adf_polar_inp_sanity_checks(adf_seed_file,cluster,adf_output_file_name,target_directory):
    #
    """Routine for the sanity checks"""
    #
    # Check the types of the files
    #
    if (type(adf_seed_file) != str):
        print('ERROR: create_adf_polar_inp called with ununderstandable adf_seed_file')
        sys.exit()
    #
    if (not isinstance(cluster,cluster_class.cluster)):
        print('ERROR: create_adf_polar_inp called without proper cluster object')
        sys.exit()
    #
    if (type(adf_output_file_name) != str):
        print('ERROR: create_adf_polar_inp called with ununderstandable adf_output_file_name')
        sys.exit()
    #
    if (type(target_directory) != str):
        print('ERROR: create_adf_polar_inp called with ununderstandable target_directory')
        sys.exit()
    #
