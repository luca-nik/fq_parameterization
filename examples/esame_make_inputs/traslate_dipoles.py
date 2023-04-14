#/home/luca/programmi/anaconda3/bin python3
# -*- coding: utf-8 -*-
import os
import sys
sys.path.insert(1, '/home/luca/programmi/des_fq_parameterization/src')
import dipoles_interface
from classes import dipoles_class
import numpy as np
#
# For each of the target DES molecule, this program:
# 1) loads the seed_dipoles;
# 2) for each dipole, and each displacement, it generates a new .dip file 
#    where only this dipole is traslated of this amount.
# 3) -optional- changes the sign to the dipoles
# 4) -optional- creates a .xyz file of the traslated dipole system.
#
molecules = ['cl', 'choline', 'ethylenglycol']
debug = False
#
# Setup the displacements
#
displacements = [float(i) for i in np.arange(2,9.75, 0.25)]
#
# Sign of the traslated dipoles
#
sign = '+-'
#
# Load the target dipoles
#
for molecule in molecules:
    dipoles = dipoles_class.dipoles()
    dipoles.initialize_from_dip('seed_dipoles/' + molecule +'_dipoles.dip')
    #
    # For each dipole and each displacement make a new .dip file and change sign if you want to
    #
    for dip in range(0, dipoles.n_dipoles):
        for j, displ, in enumerate(displacements):
            #
            new_dipoles = dipoles.move_dipoles(which_dipoles = dip, displacements = displ)
            #
            # Default sign = +-
            #
            new_dipoles.change_sign(signs = sign)
            #
            # Set the name, save into the right .xyz file, and clear the object new_dipoles
            #
            name = molecule + '_dip' + str(dip).zfill(2) +'_d' + '{:4.2f}'.format(displ+0.5).zfill(5)  
            new_dipoles.write_dip(name = name, directory = 'traslated_dipoles/'+ molecule + '/')
            #
            if (debug): #debug option, to see if the dipoles are correct
                new_dipoles.write_xyz(name = name, directory = 'traslated_dipoles/'+ molecule + '/xyz/')
            #
            new_dipoles = dipoles_class.dipoles()
