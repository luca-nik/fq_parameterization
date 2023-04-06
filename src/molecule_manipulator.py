from classes import molecule_class
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
    positions = []
    atomtypes = []
    for surf_index in molecule.surface_atoms: 
        #
        # Sanity  check
        #
        if (molecule.number_connections[surf_index] == 0):
            print('ERROR: ' + atom + str(surf_index) + ' ' + str(molecule.atomtypes(surf_index)) + \
                  ' has zero connectivity.')
            sys.exit()
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
            atomtypes.append('Mg')
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
            atomtypes.append('Mg')
            #
            versor2 = np.asarray(molecule.coords[surf_index] - molecule.coords[connector_indices[1]])/\
                      np.linalg.norm(molecule.coords[surf_index] - molecule.coords[connector_indices[1]])
            #
            positions.append(molecule.coords[surf_index] + versor2)
            atomtypes.append('Mg')
            #
            # Orthogonal dipoles
            #
            versor_ort = np.cross(versor1,versor2)
            versor_ort = versor_ort/np.linalg.norm(versor_ort)
            #
            positions.append(molecule.coords[surf_index] + versor_ort)
            atomtypes.append('Mg')
            positions.append(molecule.coords[surf_index] - versor_ort)
            atomtypes.append('Mg')
            #
            # In plane dipoles (bisectors)
            #
            versor3 = versor1 + versor2
            versor3 = versor3/np.linalg.norm(versor3)
            #
            positions.append(molecule.coords[surf_index] + versor3)
            atomtypes.append('Mg')
            positions.append(molecule.coords[surf_index] - versor3)
            atomtypes.append('Mg')
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
            atomtypes.append('Mg')
            positions.append(molecule.coords[surf_index] - versor3)
            atomtypes.append('Mg')
            #
        else:
            print("ERROR: more than three connecting atoms. Don't know what to do in this case.")
            sys.exit()

    dipoles = molecule_class.molecule(atomtypes,positions)
    return dipoles 

            
