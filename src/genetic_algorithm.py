from classes import molecule_class
from classes import dipoles_class
from classes import polarizable_embedding_class
from classes import nanofq_class
import genetic_algorithm_tools
import constants
import numpy as np
import sys
import os
import pygad #for the GA
import subprocess

def fitness_evaluator(computed_values,normalized_reference):
    #
    # Feature normalization (E-mu)/sigma
    #
    comp_energies = np.asarray(computed_values['energies'])
    comp_energies = (comp_energies - np.mean(comp_energies))/np.std(comp_energies)
    comp_polar = np.asarray(computed_values['polar'])
    comp_polar = (comp_polar - np.mean(comp_polar))/np.std(comp_polar)
    #
    loss = np.sqrt(np.linalg.norm(comp_energies-normalized_reference['energies'])**2 + \
                   np.linalg.norm(comp_polar-normalized_reference['polar'])**2)
    #
    return (1.0/loss)
#
#
#
def run_genetic_algorithm(nanofq,reference):
    #
    # Define the fitness function. Run the PE, and get its fitness
    #
    fitness_function = genetic_algorithm_tools.PE_run_and_fit
    #
    # genes = chi, eta, alpha, rq, rmu depending on the model you are training
    #
    genes = genetic_algorithm_tools.get_number_of_genes(genetic_algorithm_tools.initial_PE)
    #
    # Setup GA
    #
    sol_per_pop = 5
    elistism = sol_per_pop//4 #keep 25 % of the good boys
    #
    ga_instance = pygad.GA(num_generations = 2,                     \
                           num_parents_mating = 2,                  \
                           fitness_func=fitness_function,           \
                           sol_per_pop = sol_per_pop,               \
                           num_genes = genes,                       \
                           mutation_num_genes = genes-2,            \
                           random_mutation_min_val = 0.01,          \
                           random_mutation_max_val = 0.5,           \
                           gene_space = {'low': 0.1, 'high': 0.8},  \
                           save_best_solutions=True,                \
                           allow_duplicate_genes = False,           \
                           keep_elitism =  elistism,                \
                           stop_criteria = ["saturate_10"]          \
                           )
    #
    # Run GA
    #
    ga_instance.run()
    #
    # Take the best individual out of the last generation and make the optimal polarizable embedding
    #
    [solution,best_fit,index] = ga_instance.best_solution()#(pop_fitness = ga_instance.last_generation_fitness)
    #
    print(solution)
    print(index)
    print(best_fit)
    sys.exit()
    optimal_computed_values = genetic_algorithm_tools.final_values[index]
    #print(optimal_computed_values)
    #
    # Clean the directories
    #
    wdir = genetic_algorithm_tools.wdir
    optimal_directory = wdir+ 'g' + str(ga_instance.generations_completed)+'_p' + str(index)
    if os.path.exists(wdir + 'optimal'):
        subprocess.run(['rm', '-rf', wdir + 'optimal'])
    subprocess.run(['mv', optimal_directory, wdir + 'optimal'])
    #
    for i in range(0,sol_per_pop+1):
        try:
            subprocess.run(['rm', '-rf', wdir+ 'g' + str(ga_instance.generations_completed)+'_p'+str(i)])
        except:
            pass
    #
    # Intialize the optimal force field
    #
    optimal_embedding = polarizable_embedding_class.polarizable_embedding()
    optimal_embedding.force_field = genetic_algorithm_tools.initial_PE.force_field
    optimal_embedding.atomtypes = genetic_algorithm_tools.initial_PE.atomtypes.copy()
    optimal_embedding.pqeq = genetic_algorithm_tools.initial_PE.pqeq
    genetic_algorithm_tools.assign_new_parameters(solution,optimal_embedding)
    #
    # Print info
    #
    genetic_algorithm_tools.log_file.write('\n***************************************\n')
    genetic_algorithm_tools.log_file.write('\n')
    genetic_algorithm_tools.log_file.write('-----Optimal Polarizable Embedding-----\n')
    genetic_algorithm_tools.log_file.write('generation: ' + str(ga_instance.generations_completed) + ' member: ' + str(index) + '\n')
    genetic_algorithm_tools.log_file.write(' energy diff : ' + \
      str(np.linalg.norm(np.array(optimal_computed_values['energies'])-np.array(reference['energies']))) + '\n')
    genetic_algorithm_tools.log_file.write(' polar  diff : ' + \
      str(np.linalg.norm(np.array(optimal_computed_values['polar'])-np.array(reference['polar']))) + '\n')
    genetic_algorithm_tools.log_file.write(' fitness     : ' + str(best_fit) + '\n')
    optimal_embedding.print_info(file_=genetic_algorithm_tools.log_file)
    #
    #
    # Close the GA log_file
    #
    genetic_algorithm_tools.log_file.close()
    #
    return optimal_embedding

