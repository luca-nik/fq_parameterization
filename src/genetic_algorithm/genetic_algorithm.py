from classes import molecule_class
from classes import dipoles_class
from classes import polarizable_embedding_class
from classes import nanofq_class
from genetic_algorithm import ga_core
import constants
import numpy as np
import sys
import os
import pygad #for the GA
import subprocess
#
#
#
def run_genetic_algorithm(nanofq,reference):
    #
    # Define the fitness function. Run the PE, and get its fitness
    #
    fitness_function = ga_core.PE_run_and_fit
    #
    # Keep x% of the good people 
    #
    elistism = ga_core.ga_var['population_dimension']//ga_core.ga_var['elitism_percentage']
    #
    # run_GA
    #
    ga_instance = pygad.GA(num_generations = ga_core.ga_var["generations"],                   \
                           num_parents_mating = ga_core.ga_var["num_parents_mating"],                  \
                           fitness_func=fitness_function,           \
                           sol_per_pop = ga_core.ga_var["population_dimension"],               \
                           num_genes = ga_core.ga_var['num_genes'],                       \
                           mutation_num_genes = ga_core.ga_var["mutation_num_genes"],        \
                           random_mutation_min_val = ga_core.ga_var["mutation_min_val"],          \
                           random_mutation_max_val = ga_core.ga_var["mutation_max_val"],          \
                           gene_space = ga_core.ga_var['gene_space'],                 \
                           save_solutions=True,                     \
                           allow_duplicate_genes = False,           \
                           keep_elitism =  elistism                 \
                           #random_seed = 1,\
                           #stop_criteria = ["saturate_10"]          \
                           )
    #
    # Run GA
    #
    ga_instance.run()
    ga_instance.save('GA_state')
    #
    # Take the individuals of the last generation with their fitness
    #
    last_generation_solutions = ga_instance.solutions[-ga_core.ga_var['population_dimension']:]
    #
    best_solutions = select_best_solutions(ga_instance.keep_elitism,last_generation_solutions,\
                                           ga_instance.last_generation_fitness)
    #
    # Handle best solutions (directories, files, etc..)
    #
    best_solutions_management(ga_instance,best_solutions,ga_core.ga_var['population_dimension'])
    #
    # Close the GA log_file
    #
    ga_core.log_file.close()
    #
    # Print fitness
    #
    return best_solutions
#
##########################################################################################################################    
#
def select_best_solutions(elitism,solutions,fitness):
    #
    # Select a number of best solutions depending on the number of elitism
    #
    fitness_indices = np.flip(np.argsort(fitness))
    #
    sol = []
    fit = []
    pop_indx = []
    #
    # If no elitism keep only the best solution
    #
    if (elitism == 0):
        sol.append(solutions[0])
        fit.append(fitness  [0])
        pop_indx.append(0)
    else:
        for num, index in enumerate(fitness_indices):
            if(num == elitism):
                break
            else:
                sol.append(solutions[index])
                fit.append(fitness[index])
                pop_indx.append(index)
    #
    best_solutions = {'genes'      : sol,
                      'fitness'    : fit,
                      'pop_index'  : pop_indx}
    #
    return best_solutions
#
##########################################################################################################################    
#
def best_solutions_management(ga_instance,best_solutions,sol_per_pop):
    #
    # Manage the best solution folders and files
    #
    ga_core.log_file.write('\n***************************************\n')
    ga_core.log_file.write('\n')
    ga_core.log_file.write('-----Optimal Polarizable Embeddings----\n\n')
    polarizable_embeddings = []
    for i,genes in enumerate(best_solutions['genes']):
        ga_core.log_file.write(' member: ' + str(best_solutions["pop_index"][i]) + '\n')
        ff = ga_core.initial_PE.force_field
        pe = polarizable_embedding_class.polarizable_embedding()
        pe.force_field = ff
        pe.atomtypes  = ga_core.initial_PE.atomtypes.copy()
        pe.pqeq       = ga_core.initial_PE.pqeq
        ga_core.assign_new_parameters(genes,pe)
        polarizable_embeddings.append(pe)
        pe.print_info(file_ = ga_core.log_file)
        ga_core.log_file.write('\n')
    #
    # Manage the directories
    #
    wdir = ga_core.wdir
    for i in range(0,sol_per_pop+1):
        try:
            subprocess.run(['rm', '-rf', wdir + 'optimal_p' + str(i)])
        except:
            pass
    #
    for i,indx in enumerate(best_solutions['pop_index']):
        #
        if indx < ga_instance.keep_elitism:
            #
            # we have to remake the computation
            #
            ga_core.run_final_PE(ga_instance,dir_ = wdir + 'optimal_p' + str(indx), embedding = polarizable_embeddings[i], pop_index = indx)
        else:
            pop_indx_dir = wdir+ 'g' + str(ga_instance.generations_completed)+'_p' + str(indx)
            if (os.path.exists(pop_indx_dir)):
                subprocess.run(['mv', pop_indx_dir, wdir + 'optimal_p' + str(indx)])
            else:
                ga_core.run_final_PE(ga_instance,dir_ = wdir + 'optimal_p' + str(indx), \
                    embedding = polarizable_embeddings[i], pop_index = indx)
                
    #
    # Delete all the other folders
    #
    for i in range(0,sol_per_pop+1):
        try:
            subprocess.run(['rm', '-rf', wdir+ 'g' + str(ga_instance.generations_completed)+'_p'+str(i)])
        except:
            pass
    #
