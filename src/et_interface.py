from classes import molecule_class
from classes import dipoles_class
import numpy as np
import sys
import os
#
def create_EE_inp(et_seed_file = '', molecule_file = '',dipoles_file = '', which_dipoles = [],\
                  et_output_file_name = '', target_directory =  './', computation_name = ''):
    #
    """Procedure to generate a .inp file (input for eT calculations
       starting from a seed.inp file
       The seed.inp file, which is et_seed_file has everything set up except for:
       1) Name of the computation;
       2) Target QM molecule;
       3) Static dipole embedding the QM molecule; 
  
       The procedure is built in such a way that the user provides:
       1) et_seed_file.inp 
       2) molecule_file.xyz (target QM molecule xyz file)
       3) dipoles_file.dip  (static dipolar embedding.dip file)
       4) et_output_file_name.inp (saved in target_directory)
       5) target_directory aforementioned (skip if you already have the path in the name of the file)
       6) which_dipoles (a list of dipoles you will include in the QM/EE calculation, default is all)
       7) computation_name (an ID inserted in the .inp file, useful for the user)


       WARNING: all the files need to be provided wih the proper path to fetch them!
       WARNING: the .dip files contain information about the center of mass of the dipole
       """
       #
       # TODO
       #
       # 1) Sanity checks
       # 2) Ora metto solo il CM del dipolo, non tutto il dipolo
       # 3) Controllare che metta tutti i dipoli se richiesto
       # ALTRO
    #
    # TODO Sanity checks
    #
    pass
    #
    # Fetch the files
    #
    with open(et_seed_file, 'r') as seed:
        seed_lines = seed.readlines()
    #
    QMmolecule = molecule_class.molecule()
    QMmolecule.initialize_from_xyz(molecule_file)
    #
    EEdipoles = dipoles_class.dipoles()
    EEdipoles.initialize_from_dip(dipoles_file)
    #
    # Initialize the .inp file
    #
    with open(target_directory + et_output_file_name, 'w') as et_file:
        #
        for line in seed_lines:
            #
            if ('name:' in line):
                et_file.write('name: ' + computation_name + '\n')
            #
            if ('charge:' in line):
                et_file.write('charge: ' + str(QMmolecule.charge) + '\n')
            #
            if ('geometry' in line and 'end geometry' not in line):
                et_file.write('geometry \n')
                et_file.write('basis: aug-cc-pvdz \n')
                for i, sym in enumerate(QMmolecule.atomtypes):
                    et_file.write(sym.rjust(2) + '  ' + \
                                  '{:5.5f}'.format(QMmolecule.coords[i][0]).rjust(10)+ '  ' + \
                                  '{:5.5f}'.format(QMmolecule.coords[i][1]).rjust(10)+ '  ' + \
                                  '{:5.5f}'.format(QMmolecule.coords[i][2]).rjust(10) + '\n')
                #
                et_file.write('--')
                for i,dip in enumerate(which_dipoles):
                    et_file.write('\n' + 'H ' + ' [IMol=' + str(i+1).rjust(2) + ']  ' + \
                                   '{:5.5f}'.format(EEdipoles.positions[dip][0]).rjust(10)+ '  ' + \
                                   '{:5.5f}'.format(EEdipoles.positions[dip][1]).rjust(10)+ '  ' + \
                                   '{:5.5f}'.format(EEdipoles.positions[dip][2]).rjust(10)+ '  [q = 0.0] \n' )

                #
                # EE embedding
                #
            else:
                et_file.write(line)
    #





