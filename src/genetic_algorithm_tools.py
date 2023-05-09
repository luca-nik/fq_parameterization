from classes import molecule_class
from classes import dipoles_class
from classes import cluster_class
from classes import polarizable_embedding_class
from classes import nanofq_class
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


def global_variables_setup(workdir = '', reference_energies = [], dipoles_files = [], clusters_files = [], \
                           nanofq_seed = '', polarizable_embedding_seed = ''):
    global wdir
    global reference 
    global dip_files
    global clust_files
    global nanofq
    global initial_PE
    global log_file
    wdir = workdir
    reference = reference_energies.copy()
    dip_files = dipoles_files.copy()
    clust_files = cluster_files.copy()
    nanofq = nanofq_seed
    initial_PE = polarizable_embedding_seed
    log_file = open(wdir + 'GA_logfile.txt', 'a')
    
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
    polar  = []
    #
    # Cycle over the dipoles files to get the EE interaction energy
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
        # Setup the nanofq EE calculation
        #
        calc_name = new_nanofq.guess_name_from_dip()
        new_nanofq.name = target_directory +  '/' + calc_name
        #
        new_nanofq.create_ee_input(input_ = new_nanofq.name + '.mfq', computation_comment = new_nanofq.name, \
                                   which_dipoles = which_dipoles)
        #
        # Run it and get the energy
        #
        new_nanofq.run()
        #
        energy.append(new_nanofq.get_energy())
    #
    # Remove the directory
    #
    subprocess.run(['rm', '-rf', target_directory])
    #
    # Cycle over the clusters to get the polarizability
    #
    for clust_file in clust_files:
        #
        # Get the molecules from the selected cluster file ### the cluster object is a list of molecules, the xyz has on the second lin the way to identify the molecule
        #
        cluster = cluster_class.initialize_from_clust(clust_file)
        #
        new_nanofq = nanofq_class.nanofq(molecule = cluster, nanofq_path = nanofq.nanofq_path)
        #
        new_nanofq.polarizable_model = new_embedding
        #
        # Setup the nanofq polar
        #
        new_nanofq.name = target_directory +  clust_file.split('.clust')[0]
        #
        new_nanofq.create_polar_input(input_ = new_nanofq.name + '.mfq', computation_comment = new_nanofq.name)
        #
        # Run it and get the energy
        #
        new_nanofq.run()
        #
        # Get polar
        #
        polar.append(new_nanofq.get_polar(which = 'isotropic'))
    #
    # Remove the directory
    #
    sys.exit()
    subprocess.run(['rm', '-rf', target_directory])
    #
    # Evaluate fitness of the current individual
    #
    fitness = genetic_algorithm.fitness_function(energy,reference)
    #
    # Print some information
    #
    log_file.write('generation: ' + str(ga_instance.generations_completed) + ' member: ' + str(solution_idx) + \
                   ' energy diff: ' + str(np.linalg.norm(np.array(energy)-np.array(reference))) + '\n')
    new_embedding.print_info(file_=log_file)
    log_file.write('\n')
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
        genes = 2*len(initial_PE.atomtypes) - 1
    elif (initial_PE.force_field == 'fq_pqeq'):
        genes = 3*len(initial_PE.atomtypes) - 1
    elif (initial_PE.force_field == 'fqfmu'):
        genes = 3*len(initial_PE.atomtypes) -1 
    elif (initial_PE.force_field == 'fqfmu_pqeq'):
        genes = 5*len(initial_PE.atomtypes) -1 
    #
    return genes
#
def assign_new_parameters(GA_solution,polarizable_embedding):
    #
    # Assign the parameters to the variables of the polarizable embedding
    #
    number_atomtypes = len(polarizable_embedding.atomtypes)
    #
    polarizable_embedding.chi = []
    #
    # Set chi and eta. Zero is the chi of the hydrogen
    #
    for indx,chi in enumerate(GA_solution[0:number_atomtypes-1]):
        if (polarizable_embedding.atomtypes[indx] == 'H'):
            polarizable_embedding.chi.append(0.0)
            polarizable_embedding.chi.append(chi)
        else:
            polarizable_embedding.chi.append(chi)
    #
    # Manage the case the hydrogen is the last atomtype
    #
    if (len(polarizable_embedding.chi) < len(polarizable_embedding.atomtypes)):
        polarizable_embedding.chi.append(0.0)
    #
    polarizable_embedding.eta = [i for i in GA_solution[number_atomtypes-1:2*number_atomtypes-1]]
    #
    # Other polarizable emebdding models
    #
    if (polarizable_embedding.force_field == 'fq_pqeq'):
        polarizable_embedding.Rq  = [i for i in GA_solution[2*number_atomtypes-1:3*number_atomtypes-1]]
    #
    elif (polarizable_embedding.force_field == 'fqfmu'):
        polarizable_embedding.alpha  = [i for i in GA_solution[2*number_atomtypes-1:3*number_atomtypes-1]]
    #
    elif (polarizable_embedding.force_field == 'fqfmu_pqeq'):
        polarizable_embedding.alpha  = [i for i in GA_solution[2*number_atomtypes-1:3*number_atomtypes-1]]
        polarizable_embedding.Rq     = [i for i in GA_solution[3*number_atomtypes-1:4*number_atomtypes-1]]
        polarizable_embedding.Rmu    = [i for i in GA_solution[4*number_atomtypes-1:5*number_atomtypes-1]]
    #
def run_optimal_PE(optimal_embedding):
    #
    # Procedure to perform on the optimal PE
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
        new_nanofq.create_ee_input(input_ = new_nanofq.name + '.mfq', computation_comment = new_nanofq.name, \
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
    log_file.write('\n***************************************\n')
    log_file.write('-----Optimal Polarizable Embedding-----\n')
    optimal_embedding.print_info(file_=log_file)
    log_file.write('Optimal solution energy diff: ' + str(np.linalg.norm(np.array(energy)-np.array(reference))) + '\n')
    #
    # Evaluate fitness of the current individual
    #
    fitness = genetic_algorithm.fitness_function(energy,reference)
    #
    return fitness

