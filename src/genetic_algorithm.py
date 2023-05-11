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

def fitness_function(computed_values,reference):
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
def run_genetic_algorithm(nanofq,reference):
    #
    # Define the fitness function
    #
    fitness_function = genetic_algorithm_tools.PE_run_and_fit
    #
    # genes = chi, eta, alpha, rq, rmu depending on the model you are training
    #
    genes = genetic_algorithm_tools.get_number_of_genes(genetic_algorithm_tools.initial_PE)
    #
    #original_stdout = sys.stdout
    #warnings =  open('GA_logfile.txt', 'w')
    #sys.stdout = warnings
    ga_instance = pygad.GA(num_generations = 1,              
                           num_parents_mating = 5,            
                           fitness_func=fitness_function,   
                           sol_per_pop = 10,                 
                           num_genes = genes,               
                           mutation_num_genes = genes-2,     
                           random_mutation_min_val = 0.01,  
                           random_mutation_max_val = 0.5,   
                           gene_space = {'low': 0.1, 'high': 1}
                           )
    #sys.stdout = original_stdout
    #warnings.close()
    #
    # Select best individual and make the optimal polarizable embedding
    #
    solution = ga_instance.best_solution()[0]
    #
    # Intialize the optimal force field and make final run
    #
    optimal_embedding = polarizable_embedding_class.polarizable_embedding()
    optimal_embedding.force_field = genetic_algorithm_tools.initial_PE.force_field
    optimal_embedding.atomtypes = genetic_algorithm_tools.initial_PE.atomtypes.copy()
    optimal_embedding.pqeq = genetic_algorithm_tools.initial_PE.pqeq
    genetic_algorithm_tools.assign_new_parameters(solution,optimal_embedding)
    #
    optimal_fitness = genetic_algorithm_tools.run_optimal_PE(optimal_embedding)
    #
    # Close the GA log_file
    #
    genetic_algorithm_tools.log_file.close()
    #
    return optimal_embedding

