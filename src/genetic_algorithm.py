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
    comp_energies = np.asarray(computed_values['energies'])
    comp_polar = np.asarray(computed_values['polar'])
    #
    # Feature normalization (E-mu)/sigma
    #
    if genetic_algorithm_tools.normalization_method == 'gaussianize':
        #
        comp_energies = (comp_energies - np.mean(comp_energies))/np.std(comp_energies)
        #
        comp_polar = (comp_polar - np.mean(comp_polar))/np.std(comp_polar)
        #
        # Compute Loss sqrt(sum(E-Eref)**2/N_Eref + sum(alpha-alpha_ref)**2/N_alpha_ref
        #
        loss = np.sqrt(np.sum(np.square(comp_energies-normalized_reference['energies']))/len(comp_energies) + \
                       np.sum(np.square(comp_polar-normalized_reference['polar']))/len(comp_polar))
    #
    # Feature normalization E_i = Ecomp_i*(1/Eref_i)
    #
    else:
        #
        comp_energies = comp_energies*normalized_reference['energies']
        #
        comp_polar = np.trace(comp_polar*normalized_reference['polar'])
        #
        # Compute Loss sqrt(sum(E-1)**2/N_Eref + sum(alpha-1)**2/N_alpha_ref
        #
        loss = np.sqrt(np.sum(np.square(comp_energies - 1.0))/len(comp_energies) + \
                       np.sum(np.square(comp_polar    - 1.0))/len(comp_polar))
    
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
    sol_per_pop = 100
    elistism = sol_per_pop//4 #keep 25 % of the good boys
    #
    ga_instance = pygad.GA(num_generations = 100,                   \
                           num_parents_mating = 4,                  \
                           fitness_func=fitness_function,           \
                           sol_per_pop = sol_per_pop,               \
                           num_genes = genes,                       \
                           mutation_num_genes = genes,              \
                           random_mutation_min_val = 0.01,          \
                           random_mutation_max_val = 0.5,           \
                           gene_space = {'low': 0.1, 'high': 0.8},  \
                           save_solutions=True,                     \
                           allow_duplicate_genes = False,           \
                           keep_elitism =  elistism                 \
                           #stop_criteria = ["saturate_10"]          \
                           )
    #
    # Run GA
    #
    ga_instance.run()
    #
    # Take the individuals of the last generation with their fitness
    #
    last_generation_solutions = ga_instance.solutions[-sol_per_pop:]
    #
    best_solutions = select_best_solutions(ga_instance.keep_elitism,last_generation_solutions,\
                                           ga_instance.last_generation_fitness)
    #
    # Make a list of the optimal polarizable embeddings
    #
    genetic_algorithm_tools.log_file.write('\n***************************************\n')
    genetic_algorithm_tools.log_file.write('\n')
    genetic_algorithm_tools.log_file.write('-----Optimal Polarizable Embeddings----\n\n')
    polarizable_embeddings = []
    for i,genes in enumerate(best_solutions['genes']):
        genetic_algorithm_tools.log_file.write(' member: ' + str(best_solutions["pop_index"][i]) + '\n')
        pe = polarizable_embedding_class.polarizable_embedding()
        pe.forcefield = genetic_algorithm_tools.initial_PE.force_field
        pe.atomtypes  = genetic_algorithm_tools.initial_PE.atomtypes.copy()
        pe.pqeq       = genetic_algorithm_tools.initial_PE.pqeq
        genetic_algorithm_tools.assign_new_parameters(genes,pe)
        polarizable_embeddings.append(pe)
        pe.print_info(file_ = genetic_algorithm_tools.log_file)
        genetic_algorithm_tools.log_file.write('\n')
    #
    # Manage the directories
    #
    wdir = genetic_algorithm_tools.wdir
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
            genetic_algorithm_tools.run_single_PE(ga_instance,dir_ = wdir + 'optimal_p' + str(indx), embedding = polarizable_embeddings[i], pop_index = indx)
        else:
            pop_indx_dir = wdir+ 'g' + str(ga_instance.generations_completed)+'_p' + str(indx)
            subprocess.run(['mv', pop_indx_dir, wdir + 'optimal_p' + str(indx)])
    #
    # Delete all the other folders
    #
    for i in range(0,sol_per_pop+1):
        try:
            subprocess.run(['rm', '-rf', wdir+ 'g' + str(ga_instance.generations_completed)+'_p'+str(i)])
        except:
            pass
    #
    # Close the GA log_file
    #
    genetic_algorithm_tools.log_file.close()
    #
    return best_solutions
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
def best_solution_management(embeddings,generations):
    #
    # Here we get the generation and index of the optimal embeddings.
    # In this way we keep only the desired folders and make the computation only for the parameters not obtained in the last genration
    #
    lines = genetic_algorithm_tools.log_file.read_lines()
    #
    #for embedding in embeddings:
    #    for line in lines:
            









