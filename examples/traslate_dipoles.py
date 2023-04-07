#!/home/luca/programmi/anaconda3/bin python3
# -*- coding: utf-8 -*-
import os
import sys
sys.path.insert(1, '/home/luca/programmi/des_fq_parameterization/src')
import dipoles_interface
from classes import dipoles_class
#
# Load the target dipoles
#
dipoles = dipoles_class.dipoles()
dipoles.initialize_from_dip('test_dipoles/cl_dipoles.dip')
#
# Traslate dipoles
#
new_dipoles = dipoles.move_dipoles(which_dipoles = [0,1], displacements = [1.0, 2.0])
new_dipoles.write_xyz(name = 'cl_dipoles_traslated', directory = 'test_dipoles/')
