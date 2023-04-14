#!/home/luca/programmi/anaconda3/bin python3
# -*- coding: utf-8 -*-
import os
import sys
sys.path.insert(1, '/home/luca/programmi/des_fq_parameterization/src')
import et_interface
from classes import dipoles_class
import numpy as np
#
# For each of the target DES molecule, each generated dipole, and each distance of this dipole
# from the QM DES molecule, this program creates an eT.inp input file.
# A eT input seed file is needed. In this file you need to leave empty or with xxx:
# 1) the name of the computation;
# 2) the space for the QM and EE system;
# 3) the charge of the QM molecule. 
#    This is why you need CHARGE: keyword at the second line of the .xyz of the DES.
#
molecules = ['cl', 'choline', 'ethylenglycol']
molecules_directory = 'initial_geometries/' 
#
for molecule in molecules:
    #
    # Fetch the correct .dip files
    #
    dipoles = dipoles_class.dipoles()
    dipoles_directory = 'changed_sign_traslated_dipoles/' + molecule + '/'
    #
    files_ = [i for i in os.listdir(dipoles_directory) if i.endswith('.dip')]
    files_.sort()
    #
    # For each dipole file make the right eT input
    #
    for i,file_ in enumerate(files_):
        #
        # Get the number of the traslated dipole (python indexing)
        #
        dipole = file_.split('_')[1]
        dipole = int(dipole[-2:])
        #
        # Initialize the dipoles
        #
        dipoles.initialize_from_dip(dipoles_directory + file_)
        name = file_.split('.dip')[0]
        print(name + '.inp' + ' in the making...')
        #
        # Set up a nice name for the computation
        #
        values = name.split('_')
        dipole_number = str(int(values[1].split('dip')[-1]))
        distance = str(float(values[2].split('d')[-1]))
        #
        computation_name = values[0] + ' with embedding dipole number ' + dipole_number + \
                           ' at ' + distance + ' Ang. and sign ' + dipoles.signs[dipole]
        #
        # Create eT file (selecting the correct dipole to insert in the computation)
        #
        et_interface.create_EE_inp(et_seed_file = 'seed.inp',                                 \
                                   QMmolecule_file = molecules_directory + molecule + '.xyz', \
                                   EEdipoles_file = dipoles_directory + file_,                \
                                   which_dipoles = [dipole],                                  \
                                   et_output_file_name = name + '.inp',                       \
                                   target_directory = 'eT_inputs/' + molecule + '/',          \
                                   computation_name = computation_name)
#
