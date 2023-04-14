#!/home/luca/programmi/anaconda3/bin python3
# -*- coding: utf-8 -*-
import os
import sys
sys.path.insert(1, '/home/luca/programmi/des_fq_parameterization/src')
import dipoles_interface
from classes import molecule_class
from classes import dipoles_class
#
# For each of the target DES molecule, this program:
# 1) loads the target DES molecule;
# 2) gets its connectivity information and the atoms interacting with embedding dipoles (interface_atoms);
# 3) generates the embedding dipoles around the interface atoms;
# 4) saves the .xyz and .dip files in the folder seed_dipoles
#
molecules = ['cl', 'choline', 'ethylenglycol']
#
for molecule in molecules:
    #
    # Load the target molecule
    #
    des = molecule_class.molecule()
    des.initialize_from_xyz('initial_geometries/'+ molecule + '.xyz')
    #
    # Get the connectivity information and identify the atoms at the interface with the dipoles
    #
    des.get_connectivity(print_info    = False, bond_threshold = 1.6)
    des.get_interface_atoms(print_info = False)
    #
    # Generate the dipoles positions and save into a xyz file
    #
    dipoles = dipoles_interface.position_the_dipoles(des)
    dipoles.write_xyz(name = 'seed_dipoles/' + molecule + '_dipoles', directory = './')
    dipoles.write_dip(name = 'seed_dipoles/' + molecule + '_dipoles', directory = './')
#
