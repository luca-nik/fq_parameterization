import time
from classes import molecule_class
from classes import dipoles_class
from classes import cluster_class
from classes import polarizable_embedding_class
from classes import nanofq_class
from genetic_algorithm import genetic_algorithm
import constants
import numpy as np
import sys
import os
import pygad #for the GA
import subprocess
#
# In this subroutine you have all the core methods to run a genetic algorithm: 
# 1) global_variables_setup : set the global variables which will be passed among the different parts of the code
# 2) PE_run_and_fit         : Method to run the nanofq calcuations and handle the files
# 3) 


def global_variables_setup(workdir = '', reference_dictionary = {}, dipoles_files = [], clusters_files = [], \
                           nanofq_seed = '', polarizable_embedding_seed = '', normalization = 'to_one',\
                           input_file = ''):
    #
    """This method is required to make some global variables
    """
    #
    #
    #
    global wdir
    global reference 
    global normalized_reference 
    global dip_files
    global clust_files
    global nanofq
    global initial_PE
    global log_file
    global normalization_method 
    global ga_var
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
    if normalization == 'to_one':
        normalization_method = 'to_one'
    else:
        normalization_method = 'gaussianize'
    #
    # Feature normalization (E-mu)/sigma
    #
    if normalization_method == 'gaussianize':
        ref_energies = np.asarray(reference['energies'])
        ref_energies = (ref_energies - np.mean(ref_energies))/np.std(ref_energies)
        #
        ref_polar = np.asarray(reference['polar'])
        ref_polar = (ref_polar - np.mean(ref_polar))/np.std(ref_polar)
        #
        normalized_reference['energies'] = ref_energies.copy()
        normalized_reference['polar'] = ref_polar.copy()
    #
    # Feature normalization E = 1
    #
    else:
        ref_energies = np.asarray(reference['energies'])
        one_over_energies = 1.0/ref_energies
        #
        ref_polar = np.asarray([np.diag(np.asarray(j)) for j in reference['polar']]).flatten()
        one_over_polar = 1.0/ref_polar
        #
        normalized_reference['energies'] = one_over_energies.copy()
        normalized_reference['polar'] = one_over_polar.copy()
        #
        # Save the flattened reference polar
        #
        reference['polar'] = ref_polar
    #
    # Genetic algorithm variables setup
    #
    genes = get_number_of_genes(initial_PE)
    #
    # Setup GA variables from the input
    #
    ga_var = read_ga_variables(input_file)
    ga_var['num_genes'] = genes
    #
    #
    # Check the consistency of the genes and set up the gene_space
    #
    ga_var = setup_gene_space(ga_var)
#
##########################################################################################################################    
#
def PE_run_and_fit(ga_instance,solution,solution_idx):
    #
    """This is the main procedure done during the GA and it is done for each idividual
       1) Create a polarizable embedding with the genes of the individual
       2) Create a directory with the name of the current generation (g) and individual (p)
       3) Generate the directory for energy calculations, and make the dipoles calculations from the dip_files provided
       4) Get the energies
       5) Generate folder for polar calculatons and compute them
       6) Get polars
       7) Evaluate fitness of this individual
    """
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
        polar.append(new_nanofq.get_polar(which = 'tensor'))
    #
    # Remove the directory
    #
    if (ga_instance.generations_completed < ga_instance.num_generations):
        subprocess.run(['rm', '-rf', target_directory])
    #
    # Save the computed values in the proper way
    #
    if normalization_method == 'gaussianize':
        print('ERROR: gaussianize va sistemato nel caso which is not tensor')
        sys.exit()
    else:
        flattened_polar_diag = np.asarray([np.diag(np.asarray(j)) for j in polar]).flatten()
    #
    computed_values = {'energies': energy,
                       'polar'   : flattened_polar_diag}
    #
    # Evaluate fitness of the current individual
    #
    fitness = fitness_evaluator(computed_values,normalized_reference)
    #
    # Compute the errors with respect to the reference
    #
    energy_error = np.linalg.norm(np.array(energy)-np.array(reference['energies']))
    polar_error  = np.linalg.norm(flattened_polar_diag-reference['polar'])  
    #
    # Print some information
    #
    log_file.write('generation: ' + str(ga_instance.generations_completed) + ' member: ' + str(solution_idx) + '\n' + \
                   ' energy diff : ' + str(energy_error) + '\n')
    log_file.write(' polar  diff : ' + str(polar_error) + '\n')
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
                           ' energy diff : ' + str(energy_error) + '\n')
            info_.write(' polar  diff : ' + str(polar_error) + '\n')
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
##########################################################################################################################    
#    
def run_final_PE(ga_instance, dir_ = './', embedding = [], pop_index = 0):
    #
    # Method analogous to PE_run_and_fit, but it is needed for the last generation (handles directories in a different way)
    # Moreover it generates a new file in the final folder with the correct information
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
        polar.append(new_nanofq.get_polar(which = 'tensor'))
    #
    # Save the computed values in the proper way
    #
    if normalization_method == 'gaussianize':
        print('ERROR: gaussianize va sistemato nel caso which is not tensor')
        sys.exit()
    else:
        flattened_polar_diag = np.asarray([np.diag(np.asarray(j)) for j in polar]).flatten()
    #
    computed_values = {'energies': energy,
                       'polar'   : flattened_polar_diag}
    #
    # Evaluate fitness of the current individual
    #
    fitness = fitness_evaluator(computed_values,normalized_reference)
    #
    # Compute the errors with respect to the reference
    #
    energy_error = np.linalg.norm(np.array(energy)-np.array(reference['energies']))
    polar_error  = np.linalg.norm(flattened_polar_diag-reference['polar'])  
    #
    # Print some information
    #
    with open(target_directory + '/infofile.txt', 'w') as info_:
        info_.write('generation: ' + str(ga_instance.generations_completed) + ' member: ' + str(pop_index) + '\n' + \
                       ' energy diff : ' + str(energy_error) + '\n')
        info_.write(' polar  diff : ' + str(polar_error) + '\n')
        info_.write(' fitness     : ' + str(fitness) + '\n')
        embedding.print_info(file_=info_)
        #
        # Flush the values in the file
        #
        info_.flush()
        info_.write('\n')
    #
    # No need to retun fitness
    #
#
##########################################################################################################################    
#
def fitness_evaluator(computed_values,normalized_reference):
    #
    """Procedure to get the actual fitness value
    """
    #
    comp_energies = np.asarray(computed_values['energies'])
    comp_polar = np.asarray(computed_values['polar'])
    #
    # Feature normalization (E-mu)/sigma
    #
    if  normalization_method == 'gaussianize':
        #
        comp_energies = (comp_energies - np.mean(comp_energies))/np.std(comp_energies)
        #
        comp_polar = (comp_polar - np.mean(comp_polar))/np.std(comp_polar)
        #
        # Compute Loss sqrt(sum(E-Eref)**2/N_Eref + sum(alpha-alpha_ref)**2/N_alpha_ref
        #
        energy_error_squared = np.sum(np.square(comp_energies-normalized_reference['energies']))/len(comp_energies) 
        polar_error_squared  = np.sum(np.square(comp_polar-normalized_reference['polar']))/len(comp_polar)
        #
    #
    # Feature normalization E_i = Ecomp_i*(1/Eref_i)
    #
    else:
        #
        comp_energies = comp_energies*normalized_reference['energies']
        #
        comp_polar = computed_values['polar']*normalized_reference['polar']
        #
        # Compute Loss sqrt(sum(E-1)**2/N_Eref + sum(alpha-1)**2/N_alpha_ref
        #
        energy_error_squared = np.sum(np.square(comp_energies - 1.0))/len(comp_energies)
        polar_error_squared  = np.sum(np.square(comp_polar    - 1.0))/len(comp_polar)
        #
    #
    loss = np.sqrt(energy_error_squared + polar_error_squared)
    #
    return (1.0/loss)
#
##########################################################################################################################    
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
##########################################################################################################################    
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
##########################################################################################################################    
#
def assign_new_parameters(GA_solution,polarizable_embedding):
    #
    # Assign the parameters to the variables of the newly generated polarizable embedding
    #
    number_atomtypes = len(polarizable_embedding.atomtypes)
    #
    polarizable_embedding.chi = []
    #
    # Set chi and eta. Zero is the chi of the hydrogen
    #
    for indx,chi in enumerate(GA_solution[0:number_atomtypes-1]):
        if (polarizable_embedding.atomtypes[indx] == ga_var['atom_constrained']):
            polarizable_embedding.chi.append(ga_var['chi_constrained'])
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
##########################################################################################################################    
#
def read_ga_variables(input_file):
    #
    #
    #
    with open(input_file, 'r') as inp_:
        lines = inp_.readlines()
    #
    ga_var = {'generations'          : 0   ,
              'num_parents_mating'   : 0   ,
              'population_dimension' : 0   ,
              'elitism_percentage'   : 5   ,
              'num_genes'            : 0   ,
              'mutation_num_genes'   : 0   , 
              'mutation_min_val'     : 0.0 ,
              'mutation_max_val'     : 0.0 ,
              'atom_constrained'     : ''  ,
              'chi_constrained'      : 0.0 ,
              'at_types_input'       : []  ,
              'chis_range'           : []  ,
              'etas_range'           : []  ,
              'alphas_range'         : []  ,  
              'gene_space'           : []
              }
    #
    #
    for line in lines:
        values = line.strip().split('=')
        if(len(values) > 1):
            #
            if (values[0] == 'generations'          or \
                values[0] == 'num_parents_mating'   or \
                values[0] == 'population_dimension' or \
                values[0] == 'mutation_num_genes'      ):
                ga_var[values[0]] = int(values[1])
            #
            elif(values[0] == 'mutation_min_val' or \
                 values[0] == 'mutation_max_val' or \
                 values[0] == 'elitism_percentage'):
                ga_var[values[0]] = float(values[1])
            #
            elif(values[0] == 'atom'):
                parameters = values[1].split(',')
                at_type = parameters[0].strip().upper()
                ga_var['at_types_input'].append(at_type)
                #
                # chi
                #
                chis    = parameters[1]
                if (len(chis.split('-')) > 1):
                    ga_var['chis_range'].append([float(i) for i in chis.split('-')[0:2]])
                else:
                    if (ga_var['atom_constrained'] != ''):
                        print('ERROR: ' + input_file + ' called with more than one constraint.')
                        sys.exit()
                    #
                    ga_var['chis_range'].append([float(chis)])
                    ga_var['atom_constrained'] = at_type
                    ga_var['chi_constrained'] = float(chis)
                #
                # etas
                #
                etas    = parameters[2]
                if (len(etas.split('-')) > 1):
                    ga_var['etas_range'].append([float(i) for i in etas.split('-')[0:2]])
                else:
                    print('ERROR: ' + input_file + ' provided with wrong eta range for ' + at_type)
                    sys.exit()
                #
                # alphas
                #
                if (initial_PE.force_field == 'fqfmu' or \
                    initial_PE.force_field == 'fqfmu_pqeq'):
                   if len(parameters) < 4:
                       print('ERROR: ' + input_file + ' provided for a fqfmu calculation without correct number of parameters')
                       sys.exit()
                   else:
                       alphas = parameters[3]
                       if (len(alphas.split('-')) > 1):
                           ga_var['alphas_range'].append([float(i) for i in alphas.split('-')[0:2]])
                       else:
                           print('ERROR: ' + input_file + ' provided with wrong alpha range for ' + at_type)
                           sys.exit()

    #
    # You need a zero of the electronegativity
    #
    if ga_var['atom_constrained'] == '':
        print('ERROR: '+ input_file + ' provided without a constraint on the chi')
        sys.exit()
    #
    return ga_var
#
##########################################################################################################################    
#
def setup_gene_space(ga_var):
   #
   """This procedure is needed to set uo the gene space.
      In particular to limit the range of variability of each gene as set up from the input
      and to put in the right order (the order of the atomtypes of the polarizable embedding)
   """
   #
   at     = ga_var['at_types_input'].copy()
   chis   = ga_var['chis_range'].copy()
   etas   = ga_var['etas_range'].copy()
   alphas = ga_var['alphas_range'].copy()
   #
   ga_var['at_types_input'] = []
   ga_var['chis_range']     = []
   ga_var['etas_range']     = []
   ga_var['alphas_range']   = []
   #
   # Rearrange the atomtypes according to those of the polarizable emebedding
   # In this way we avoid possible issues of wrong atomtypes ordering from the input file
   #
   for i, at_type in enumerate(initial_PE.atomtypes):
       for j, at_type_ in enumerate(at):
           if at_type_.lower() == at_type.lower():
               #
               ga_var['at_types_input'].append(at[j].upper())
               ga_var['chis_range'].append(chis[j])          
               ga_var['etas_range'].append(etas[j])
               #
               if (initial_PE.force_field == 'fqfmu' or \
                   initial_PE.force_field == 'fqfmu_pqeq'):
                   ga_var['alphas_range'].append(alphas[j])
   #
   gene_space = []
   #
   # Chi gene space setup
   #
   for i, chis in enumerate(ga_var['chis_range']):
       #
       if ga_var['atom_constrained'] == ga_var['at_types_input'][i]:
           pass
       else:
           low_chi  = chis[0] 
           high_chi = chis[1] 
           gene_space.append({'low' : low_chi, 'high' : high_chi})
   #
   # Chi gene space setup
   #
   for i, etas in enumerate(ga_var['etas_range']):
       #
       low_eta  = etas[0] 
       high_eta = etas[1] 
       gene_space.append({'low' : low_eta, 'high' : high_eta})
   #
   # alpha gene space setup
   #
   if (initial_PE.force_field == 'fqfmu' or \
       initial_PE.force_field == 'fqfmu_pqeq'):
       for i, alphas in enumerate(ga_var['alphas_range']):
           #
           low_alpha  = alphas[0] 
           high_alpha = alphas[1] 
           gene_space.append({'low' : low_alpha, 'high' : high_alpha})
   #
   ga_var['gene_space'] = gene_space
   # 
   return ga_var
