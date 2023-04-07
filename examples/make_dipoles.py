#!/home/luca/programmi/anaconda3/bin python3
# -*- coding: utf-8 -*-
import os
import sys
sys.path.insert(1, '/home/luca/programmi/des_fq_parameterization/src')
import dipoles_interface
from classes import molecule_class
from classes import dipoles_class
#
# Load the target molecule
#
ethylenglycol = molecule_class.molecule()
ethylenglycol.initialize_from_xyz('geometries/ethylenglycol.xyz')
#
# Get the connectivity information and identify the atoms at the interface with the DES
#
ethylenglycol.get_connectivity(print_info    = True, bond_threshold = 1.6)
ethylenglycol.get_interface_atoms(print_info = True)
#
# Generate the dipoles positions and save into a xyz file
#
dipoles = dipoles_interface.position_the_dipoles(ethylenglycol)
dipoles.write_xyz(name = 'test_dipoles/ethylenglycol_dipoles', directory = './')
#
choline = molecule_class.molecule()
choline.initialize_from_xyz('geometries/choline.xyz')
#
# Get the connectivity information and identify the atoms at the interface with the DES
#
choline.get_connectivity(print_info    = True, bond_threshold = 1.6)
choline.get_interface_atoms(print_info = True)
#
# Generate the dipoles positions and save into a xyz file
#
dipoles = dipoles_interface.position_the_dipoles(choline)
dipoles.write_xyz(name = 'test_dipoles/choline_dipoles', directory = './')
#
cl = molecule_class.molecule()
cl.initialize_from_xyz('geometries/cl.xyz')
#
# Get the connectivity information and identify the atoms at the interface with the DES
#
cl.get_connectivity(print_info    = True, bond_threshold = 1.6)
cl.get_interface_atoms(print_info = True)
#
# Generate the dipoles positions and save into a xyz file
#
dipoles = dipoles_interface.position_the_dipoles(cl)
dipoles.write_xyz(name = 'test_dipoles/cl_dipoles', directory = './')
dipoles.write_dip(name = 'test_dipoles/cl_dipoles', directory = './')
#

