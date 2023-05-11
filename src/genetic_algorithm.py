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

def fitness_evaluator(computed_values,reference):
    #
    # Feature normalization (E-mu)/sigma
    #
    comp_energies = np.asarray(computed_values['energies'])
    comp_energies = (comp_energies - np.mean(comp_energies))/np.std(comp_energies)
    comp_polar = np.asarray(computed_values['polar'])
    comp_polar = (comp_polar - np.mean(comp_polar))/np.std(comp_polar)
    #
    loss = np.sqrt(np.linalg.norm(comp_energies-reference['energies'])**2 + \
                   np.linalg.norm(comp_polar-reference['polar'])**2)
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
    sol_per_pop = 3
    ga_instance = pygad.GA(num_generations = 1,                     \
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
                           stop_criteria = ["saturate_10"]          \
                           )
    #
    # Run GA
    #
    ga_instance.run()
    #
    # Select best individual and make the optimal polarizable embedding
    #
    solutions = ga_instance.best_solutions
    fit = ga_instance.best_solutions_fitness
    index = fit.index(np.max(fit))
    #
    solution = solutions[index]
    #
    # Clean
    #
    wdir = genetic_algorithm_tools.wdir
    optimal_directory = wdir+ 'g' + str(ga_instance.generations_completed)+'_p' + str(index+1)
    subprocess.run(['mv', optimal_directory, 'optimal'])
    for i in range(1,sol_per_pop):
        try:
            subprocess.run(['rm', '-rf', wdir+ 'g' + str(ga_instance.generations_completed)+'_p'+str(i)])
        except:
            pass
    #
    # Intialize the optimal force field and make final run
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
    genetic_algorithm_tools.log_file.write('-----Optimal Polarizable Embedding-----\n')
    genetic_algorithm_tools.log_file.write('generation: ' + str(ga_instance.generations_completed) + ' member: ' + str(index) + '\n')
    optimal_embedding.print_info(file_=genetic_algorithm_tools.log_file)
    #
    #optimal_fitness = genetic_algorithm_tools.run_optimal_PE(optimal_embedding)
    #
    # Close the GA log_file
    #
    genetic_algorithm_tools.log_file.close()
    #
    return optimal_embedding

