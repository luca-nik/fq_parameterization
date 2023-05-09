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

def fitness_function(energy,reference):
    return (1.0/np.linalg.norm(np.array(energy)-np.array(reference)))
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
    ga_instance = pygad.GA(num_generations=50,              
                           num_parents_mating=5,            
                           fitness_func=fitness_function,   
                           sol_per_pop=100,                 
                           num_genes = genes,               
                           mutation_num_genes = genes-2,     
                           random_mutation_min_val = 0.01,  
                           random_mutation_max_val = 0.5,   
                           gene_space = {'low': 0, 'high': 1}
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

