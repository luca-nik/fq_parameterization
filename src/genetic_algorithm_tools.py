import time
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


def global_variables_setup(workdir = '', reference_dictionary = {}, dipoles_files = [], clusters_files = [], \
                           nanofq_seed = '', polarizable_embedding_seed = ''):
    global wdir
    global reference 
    global normalized_reference 
    global dip_files
    global clust_files
    global nanofq
    global initial_PE
    global log_file
    #
    wdir = workdir
    reference = reference_dictionary.copy()
    normalized_reference = reference_dictionary.copy()
    dip_files = dipoles_files.copy()
    clust_files = clusters_files.copy()
    nanofq = nanofq_seed
    initial_PE = polarizable_embedding_seed
    log_file = open(wdir + 'GA_logfile.txt', 'w+')
    #
    # Feature normalization (E-mu)/sigma
    #
    ref_energies = np.asarray(reference['energies'])
    ref_energies = (ref_energies - np.mean(ref_energies))/np.std(ref_energies)
    ref_polar = np.asarray(reference['polar'])
    ref_polar = (ref_polar - np.mean(ref_polar))/np.std(ref_polar)
    normalized_reference['energies'] = ref_energies.copy()
    normalized_reference['polar'] = ref_polar.copy()
    
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
    #
    try:
        os.mkdir(target_directory)
    except FileExistsError:
        subprocess.run(['rm', '-rf', target_directory])
        os.mkdir(target_directory)
    #
    os.mkdir(target_directory + '/energies')
    print('g' + str(ga_instance.generations_completed)+'_p' + str(solution_idx))
    print(ga_instance.population)
    print(solution)
    print('')
    #
    energy = []
    polar  = []
    #
    # Cycle over the dipoles files to get the EE interaction energy
    #
    for dip_file in dip_files[0:4]:
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
        new_nanofq.name = target_directory +  '/energies/' + calc_name
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
    # Restore the directory
    #
    if (ga_instance.generations_completed < ga_instance.num_generations):
        subprocess.run(['rm', '-rf', target_directory + '/energies'])
    #
    os.mkdir(target_directory + '/polar')
    #
    # Cycle over the clusters to get the polarizability
    #
    for clust_file in clust_files:
        #
        #
        #
        cluster = cluster_class.cluster()
        cluster.initialize_from_clust(clust_file)
        #
        new_nanofq = nanofq_class.nanofq()
        new_nanofq.molecule = cluster
        new_nanofq.nanofq_path = nanofq.nanofq_path
        #
        new_nanofq.polarizable_model = new_embedding
        #
        # Setup the nanofq polar
        #
        clust_name = clust_file.split('/')[-1]
        new_nanofq.name = target_directory + '/polar/' + clust_name.split('.clust')[0]
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
    if (ga_instance.generations_completed < ga_instance.num_generations):
        subprocess.run(['rm', '-rf', target_directory])
        
    #
    # Evaluate fitness of the current individual
    #
    computed_values = {'energies': energy,
                       'polar'   : polar}
    #
    fitness = genetic_algorithm.fitness_evaluator(computed_values,normalized_reference)
    #
    # Print some information
    #
    log_file.write('generation: ' + str(ga_instance.generations_completed) + ' member: ' + str(solution_idx) + '\n' + \
                   ' energy diff : ' + str(np.linalg.norm(np.array(energy)-np.array(reference['energies']))) + '\n')
    log_file.write(' polar  diff : ' + str(np.linalg.norm(np.array(polar)-np.array(reference['polar']))) + '\n')
    log_file.write(' fitness     : ' + str(fitness) + '\n')
    new_embedding.print_info(file_=log_file)
    #
    # Flush the values in the file
    #
    log_file.flush()
    log_file.write('\n')
    #
    # Print PE specific information at the last cycle in the folder
    #
    if (ga_instance.generations_completed == ga_instance.num_generations):
        with open(target_directory + '/infofile.txt', 'w') as info_:
            info_.write('generation: ' + str(ga_instance.generations_completed) + ' member: ' + str(solution_idx) + '\n' + \
                           ' energy diff : ' + str(np.linalg.norm(np.array(energy)-np.array(reference['energies']))) + '\n')
            info_.write(' polar  diff : ' + str(np.linalg.norm(np.array(polar)-np.array(reference['polar']))) + '\n')
            info_.write(' fitness     : ' + str(fitness) + '\n')
            new_embedding.print_info(file_=info_)
            #
            # Flush the values in the file
            #
            info_.flush()
            info_.write('\n')
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
#
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
#
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
#    
#    
def run_single_PE(ga_instance, dir_ = './', embedding = [], pop_index = 0):
    #
    # 
    target_directory = dir_
    #
    try:
        os.mkdir(target_directory)
    except FileExistsError:
        subprocess.run(['rm', '-rf', target_directory])
        os.mkdir(target_directory)
    #
    os.mkdir(target_directory + '/energies')
    #
    energy = []
    polar  = []
    #
    # Cycle over the dipoles files to get the EE interaction energy
    #
    for dip_file in dip_files[0:4]:
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
        new_nanofq.polarizable_model = embedding
        #
        # Setup the nanofq EE calculation
        #
        calc_name = new_nanofq.guess_name_from_dip()
        new_nanofq.name = target_directory +  '/energies/' + calc_name
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
    # Make polar
    #
    os.mkdir(target_directory + '/polar')
    #
    # Cycle over the clusters to get the polarizability
    #
    for clust_file in clust_files:
        #
        #
        #
        cluster = cluster_class.cluster()
        cluster.initialize_from_clust(clust_file)
        #
        new_nanofq = nanofq_class.nanofq()
        new_nanofq.molecule = cluster
        new_nanofq.nanofq_path = nanofq.nanofq_path
        #
        new_nanofq.polarizable_model = embedding
        #
        # Setup the nanofq polar
        #
        clust_name = clust_file.split('/')[-1]
        new_nanofq.name = target_directory + '/polar/' + clust_name.split('.clust')[0]
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
    # Evaluate fitness of the current individual
    #
    computed_values = {'energies': energy,
                       'polar'   : polar}
    #
    fitness = genetic_algorithm.fitness_evaluator(computed_values,normalized_reference)
    #
    # Print some information
    #
    with open(target_directory + '/infofile.txt', 'w') as info_:
        info_.write('generation: ' + str(ga_instance.generations_completed) + ' member: ' + str(pop_index) + '\n' + \
                       ' energy diff : ' + str(np.linalg.norm(np.array(energy)-np.array(reference['energies']))) + '\n')
        info_.write(' polar  diff : ' + str(np.linalg.norm(np.array(polar)-np.array(reference['polar']))) + '\n')
        info_.write(' fitness     : ' + str(fitness) + '\n')
        embedding.print_info(file_=info_)
        #
        # Flush the values in the file
        #
        info_.flush()
        info_.write('\n')
    #
#
