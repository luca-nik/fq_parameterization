from classes import molecule_class
from classes import dipoles_class
from classes import polarizable_embedding_class
from classes import nanofq_class
import nanofq_interface
import genetic_algorithm
import constants
import numpy as np
import sys
import os
import pygad #for the GA
import subprocess
#
# Necessito di un get_reference
# Necessito di una classe comoda, per definire il posto dove lavoro, i dip files su cui devo andare a lavorare
#


def global_variables_setup(workdir = '', reference_energies = [], dipoles_files = [], nanofq_seed = '', polarizable_embedding_seed = ''):
    global wdir
    global reference 
    global dip_files
    global nanofq
    global initial_PE
    global log_file
    wdir = workdir
    reference = reference_energies
    dip_files = dipoles_files
    nanofq = nanofq_seed
    initial_PE = polarizable_embedding_seed
    log_file = open(wdir + 'GA_logfile.txt', 'w')
#
#
def PE_run_and_fit(ga_instance,solution,solution_idx):
    #
    """This is the main procedure done during the GA and it is done for each idividual
       1) Create a polarizable embedding with the genes of the individual
       2) Create a directory with the name of the current generation (g) and individual (p)"""
    #
    # Initialize the new polarizable embedding of this individual
    # 
    new_embedding = polarizable_embedding_class.polarizable_embedding()
    new_embedding.force_field = initial_PE.force_field
    new_embedding.atomtypes = initial_PE.atomtypes.copy()
    new_embedding.pqeq = initial_PE.pqeq
    #
    # Assign the correct new parameters depending on the model you are parameterizing
    #
    assign_new_parameters(solution,new_embedding)
    # 
    target_directory = wdir+ 'g' + str(ga_instance.generations_completed)+'_p' + str(solution_idx)
    os.mkdir(target_directory)
    #
    energy = []
    #
    # Cycle over the dipoles files
    #
    for dip_file in dip_files:
        #
        # Initialize dipoles and get the dipole you need to place
        #
        dipoles = dipoles_class.dipoles()
        dipoles.initialize_from_dip(dip_file)
        #
        which_dipoles = get_which_dipoles_from_dip(dip_file)
        #
        # Create a nano_fq object with the selected molecule, dipoles and the common path and the selected polarizable embedding
        #
        new_nanofq = nanofq_class.nanofq(molecule = nanofq.molecule, dipoles = dipoles, nanofq_path = nanofq.nanofq_path)
        #
        new_nanofq.which_dipoles = which_dipoles.copy()
        new_nanofq.polarizable_model = new_embedding
        #
        # Setup the nanofq calculation
        #
        calc_name = new_nanofq.guess_name_from_dip()
        new_nanofq.name = target_directory +  '/' + calc_name
        #
        new_nanofq.create_input(input_ = new_nanofq.name + '.mfq', computation_comment = new_nanofq.name, \
                                which_dipoles = which_dipoles)
        #
        # Run it and get the energy
        #
        new_nanofq.run()
        #
        energy.append(new_nanofq.get_energy())
    #
    # Print some information
    #
    log_file.write('generation: ' + str(ga_instance.generations_completed) + ' member: ' + str(solution_idx) + \
                   ' energy diff: ' + str(np.linalg.norm(np.array(energy)-np.array(reference))) + '\n')
    #
    # Evaluate fitness of the current individual
    #
    fitness = genetic_algorithm.fitness_function(energy,reference)
    #
    # Remove the directory
    #
    subprocess.run(['rm', '-rf', target_directory])
    #
    return fitness
#
#
#
def get_which_dipoles_from_dip(dipfile):
    #
    # Get the dipole you need to nsert in the calculation
    #
    name = dipfile.split('.dip')[0]
    name = name.split('/')[-1]
    name = name.split('dip')[1]
    name = name.split('_')[0]
    return [int(name)]
#
def get_number_of_genes(initial_PE):
    #
    # Get the number of genes dependign on the class of the PE
    #
    if (initial_PE.force_field == 'fq'):
        genes = 2*len(initial_PE.atomtypes)
    elif (initial_PE.force_field == 'fq_pqeq'):
        genes = 3*len(initial_PE.atomtypes)
    elif (initial_PE.force_field == 'fqfmu'):
        genes = 3*len(initial_PE.atomtypes)
    elif (initial_PE.force_field == 'fqfmu_pqeq'):
        genes = 5*len(initial_PE.atomtypes)
    #
    return genes
#
def assign_new_parameters(GA_solution,polarizable_embedding):
    #
    # Assign the parameters to the variables of the polarizable embedding
    #
    if (polarizable_embedding.force_field == 'fq'):
        polarizable_embedding.chi = [i for i in GA_solution[0:2]]
        polarizable_embedding.eta = [i for i in GA_solution[2:4]]
    elif (polarizable_embedding.force_field == 'fq_pqeq'):
        polarizable_embedding.chi = [i for i in GA_solution[0:2]]
        polarizable_embedding.eta = [i for i in GA_solution[2:4]]
        polarizable_embedding.Rq  = [i for i in GA_solution[4:6]]
    elif (polarizable_embedding.force_field == 'fqfmu'):
        polarizable_embedding.chi    = [i for i in GA_solution[0:2]]
        polarizable_embedding.eta    = [i for i in GA_solution[2:4]]
        polarizable_embedding.alpha  = [i for i in GA_solution[4:6]]
    elif (polarizable_embedding.force_field == 'fqfmu_pqeq'):
        polarizable_embedding.chi    = [i for i in GA_solution[0:2]]
        polarizable_embedding.eta    = [i for i in GA_solution[2:4]]
        polarizable_embedding.alpha  = [i for i in GA_solution[4:6]]
        polarizable_embedding.Rq     = [i for i in GA_solution[6:8]]
        polarizable_embedding.Rm     = [i for i in GA_solution[8:10]]
    #
def run_optimal_PE(optimal_embedding):
    #
    target_directory = wdir+ 'optimal'
    #
    os.mkdir(target_directory)
    #
    energy = []
    #
    # Cycle over the dipoles files
    #
    for dip_file in dip_files:
        #
        # Initialize dipoles and get the dipole you need to place
        #
        dipoles = dipoles_class.dipoles()
        dipoles.initialize_from_dip(dip_file)
        #
        which_dipoles = get_which_dipoles_from_dip(dip_file)
        #
        # Create a nano_fq object with the selected molecule, dipoles and the common path and the selected polarizable embedding
        #
        new_nanofq = nanofq_class.nanofq(molecule = nanofq.molecule, dipoles = dipoles, nanofq_path = nanofq.nanofq_path)
        #
        new_nanofq.which_dipoles = which_dipoles.copy()
        new_nanofq.polarizable_model = optimal_embedding
        #
        # Setup the nanofq calculation
        #
        calc_name = new_nanofq.guess_name_from_dip()
        new_nanofq.name = target_directory +  '/' + calc_name
        #
        new_nanofq.create_input(input_ = new_nanofq.name + '.mfq', computation_comment = new_nanofq.name, \
                                which_dipoles = which_dipoles)
        #
        # Run it and get the energy
        #
        new_nanofq.run()
        #
        energy.append(new_nanofq.get_energy())
    #
    # Print some information
    #
    log_file.write('\n-----Optimal Polarizable Embedding-----\n')
    log_file.write('force_field : ' + optimal_embedding.force_field + '\n')
    infostring = ''
    for ii,i in enumerate(optimal_embedding.atomtypes):
        infostring += "'" + i + "'"
        if ii < len(optimal_embedding.atomtypes) -1:
            infostring += ' , '
    log_file.write('atomtypes   : ' + infostring + '\n')
    #
    infostring = ''
    for ii,i in enumerate(optimal_embedding.chi):
        infostring += str(i)
        if ii < len(optimal_embedding.chi) -1:
            infostring += ' , '
    log_file.write('chi         : ' + infostring + '\n')
    #
    infostring = ''
    for ii,i in enumerate(optimal_embedding.eta):
        infostring += str(i)
        if ii < len(optimal_embedding.eta) -1:
            infostring += ' , '
    log_file.write('eta         : ' + infostring + '\n')
    #
    infostring = ''
    for ii,i in enumerate(optimal_embedding.alpha):
        infostring += str(i)
        if ii < len(optimal_embedding.alpha) -1:
            infostring += ' , '
    log_file.write('alpha       : ' + infostring + '\n')
    #
    infostring = ''
    for ii,i in enumerate(optimal_embedding.Rq):
        infostring += str(i)
        if ii < len(optimal_embedding.Rq) -1:
            infostring += ' , '
    log_file.write('Rq          : ' + infostring + '\n')
    #
    infostring = ''
    for ii,i in enumerate(optimal_embedding.Rmu):
        infostring += str(i)
        if ii < len(optimal_embedding.Rmu) -1:
            infostring += ' , '
    log_file.write('Rmu         : ' + infostring + '\n')
    log_file.write('Optimal solution energy diff: ' + str(np.linalg.norm(np.array(energy)-np.array(reference))) + '\n')
    #
    # Evaluate fitness of the current individual
    #
    fitness = genetic_algorithm.fitness_function(energy,reference)
    #
    return fitness

