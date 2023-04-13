from classes import molecule_class
from classes import dipoles_class
import numpy as np
import sys
import os
#
def position_the_dipoles(molecule, print_info = False):
    #
    """Procedure to generate the location of the fixed dipoles base on: """
    """ -The connectivity of the molecule
        -Its surface atoms
        -Other stuff"""
    #
    n_dipoles = 0
    positions = []
    versors   = []
    #
    for surf_index in molecule.surface_atoms: 
        #
        # Sanity  check
        #
        if (molecule.number_connections[surf_index] == 0 and molecule.atomtypes[surf_index] != 'Cl'):
            print('ERROR: ' + atom + str(surf_index) + ' ' + str(molecule.atomtypes(surf_index)) + \
                  ' has zero connectivity.')
            sys.exit()
        #
        # If its a Cl atom, then generate the dipoles along the lone pairs directions
        #
        elif (molecule.number_connections[surf_index] == 0 and molecule.atomtypes[surf_index] == 'Cl'):
            #
            coords = molecule.coords[surf_index].copy()
            #
            A_coordinates = np.asarray([coords[0] + 1.0, coords[1] + 1.0, coords[2] + 1.0])
            B_coordinates = np.asarray([coords[0] - 1.0, coords[1] - 1.0, coords[2] + 1.0])
            C_coordinates = np.asarray([coords[0] - 1.0, coords[1] + 1.0, coords[2] - 1.0])
            D_coordinates = np.asarray([coords[0] + 1.0, coords[1] - 1.0, coords[2] - 1.0])
            #
            A_versor = np.asarray(molecule.coords[surf_index]     - A_coordinates)/\
                       np.linalg.norm(molecule.coords[surf_index] - A_coordinates)
            B_versor = np.asarray(molecule.coords[surf_index]     - B_coordinates)/\
                       np.linalg.norm(molecule.coords[surf_index] - B_coordinates)
            C_versor = np.asarray(molecule.coords[surf_index]     - C_coordinates)/\
                       np.linalg.norm(molecule.coords[surf_index] - C_coordinates)
            D_versor = np.asarray(molecule.coords[surf_index]     - D_coordinates)/\
                       np.linalg.norm(molecule.coords[surf_index] - D_coordinates)
            #
            positions.append(molecule.coords[surf_index] + A_versor)
            versors.append(A_versor)
            positions.append(molecule.coords[surf_index] + B_versor)
            versors.append(B_versor)
            positions.append(molecule.coords[surf_index] + C_versor)
            versors.append(C_versor)
            positions.append(molecule.coords[surf_index] + D_versor)
            versors.append(D_versor)
            n_dipoles += 4
            #
        #
        elif (molecule.number_connections[surf_index] == 1):
            #
            # Generate dipoles position:
            # 1) Along the line connceting it to the neighbours
            #
            connector_index = molecule.connected_to[surf_index][0]
            versor = np.asarray(molecule.coords[surf_index] - molecule.coords[connector_index])/\
                     np.linalg.norm(molecule.coords[surf_index] - molecule.coords[connector_index])
            #
            positions.append(molecule.coords[surf_index] + versor)
            versors.append(versor)
            n_dipoles += 1
            #
        elif (molecule.number_connections[surf_index] == 2):
            #
            # Generate dipoles position:
            # 1) Along the line connecting it to the neighbours
            # 2) Along the visector of the angle formed with the two neighbours
            # 3) Along the direction orthogonal to the plane identified by these atoms
            #
            connector_indices = molecule.connected_to[surf_index][:]
            #
            # Directions and in plane dipoles along connections
            #
            versor1 = np.asarray(molecule.coords[surf_index] - molecule.coords[connector_indices[0]])/\
                      np.linalg.norm(molecule.coords[surf_index] - molecule.coords[connector_indices[0]])
            #
            positions.append(molecule.coords[surf_index] + versor1)
            versors.append(versor1)
            n_dipoles += 1
            #
            versor2 = np.asarray(molecule.coords[surf_index] - molecule.coords[connector_indices[1]])/\
                      np.linalg.norm(molecule.coords[surf_index] - molecule.coords[connector_indices[1]])
            #
            positions.append(molecule.coords[surf_index] + versor2)
            versors.append(versor2)
            n_dipoles += 1
            #
            # Orthogonal dipoles
            #
            versor_ort = np.cross(versor1,versor2)
            versor_ort = versor_ort/np.linalg.norm(versor_ort)
            #
            positions.append(molecule.coords[surf_index] + versor_ort)
            versors.append(versor_ort)
            n_dipoles += 1
            #
            positions.append(molecule.coords[surf_index] - versor_ort)
            versors.append(-versor_ort)
            n_dipoles += 1
            #
            # In plane dipoles (bisectors)
            #
            versor3 = versor1 + versor2
            versor3 = versor3/np.linalg.norm(versor3)
            #
            positions.append(molecule.coords[surf_index] + versor3)
            versors.append(versor3)
            n_dipoles += 1
            #
            #positions.append(molecule.coords[surf_index] - versor3) #this is on the bisector but on the acute side
            # 
        elif (molecule.number_connections[surf_index] == 3):
            #
            # Generate dipoles position symmetrically over and below the plane identified by the atom and
            # the other atoms it is connected to. We try to use the backbone of the molecule to build the 
            # planes
            #
            connector_indices = molecule.connected_to[surf_index][:]
            connectors_not_on_surface = [i for i in connector_indices if i not in molecule.surface_atoms]
            #
            if (len(connectors_not_on_surface) != 2):
                #
                #Either 3 or less than 2 remaining neighbours: i take the first two to define the plane
                #
                connector_indices = connector_indices[0:2]
            else:
                connector_indices = connectors_not_on_surface.copy()
            #
            # Now generate the direction orthogonal to the plane
            #
            versor1 = np.asarray(molecule.coords[surf_index] - molecule.coords[connector_indices[0]])/\
                      np.linalg.norm(molecule.coords[surf_index] - molecule.coords[connector_indices[0]])
            #
            versor2 = np.asarray(molecule.coords[surf_index] - molecule.coords[connector_indices[1]])/\
                      np.linalg.norm(molecule.coords[surf_index] - molecule.coords[connector_indices[1]])
            #
            versor3 = np.cross(versor1,versor2)
            versor3 = versor3/np.linalg.norm(versor3)
            #
            positions.append(molecule.coords[surf_index] + versor3)
            versors.append(versor3)
            n_dipoles += 1
            positions.append(molecule.coords[surf_index] - versor3)
            versors.append(-versor3)
            n_dipoles += 1
            #
        else:
            print("ERROR: more than three connecting atoms. Don't know what to do in this case.")
            sys.exit()

    dipoles = dipoles_class.dipoles(n_dipoles = n_dipoles, positions = positions, directions = versors)
    return dipoles 
