import sys
#
def number_connections():
    """A function which returns a dictionary."""
    """The values indicate the threshold for the connectrivity"""
    connections = {
    "H"  :   1        ,
    "HW" :   1        ,
    "C"  :   3        ,
    "N"  :   3        ,
    "O"  :   2        ,
    "OW" :   2        ,
    "Cl" :   4        ,
    }

    return connections
#
def dipoles_distance():
    return 0.5 # in Angstrom
#
def FQ_parameters(solvent = 'wat'):
    #
    # FQ parameters taken from Ambrosetti Skoko
    #    
    if solvent == 'acn':
        atomtypes = ['C', 'N', 'H']
        #
        chi = [0.04990100, 0.16615000, 0.01069500]
        eta = [0.37625000, 0.44569000, 0.57298000]
    elif solvent == 'dio' or solvent == 'thf':
        atomtypes = ['O', 'C', 'H']
        chi = [0.49112000, 0.43733000, 0.35080000]
        eta = [0.46457000, 0.39730000, 0.58406000]
    elif solvent == 'met' or solvent == 'eth':
        atomtypes = ['C', 'O', 'H']
        chi = [7.6251e-01, 9.5000e-01, 5.8328e-01]
        eta = [1.7382e-01, 9.5000e-01, 5.2209e-01]
    elif solvent == 'wat':
        atomtypes = ['O', 'H']
        chi = [0.15780000, 0.01249200]
        eta = [0.54386000, 0.44020000]
    else:
        print('Error: I do not have FQ parameters for ' + str(solvent))
        sys.exit()
    return atomtypes, chi, eta
 

