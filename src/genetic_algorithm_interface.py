from classes import molecule_class
from classes import dipoles_class
from classes import polarizable_embedding_class
from classes import nanofq_class
import nanofq_interface
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


def global_variables_setup(workdir = '', reference_energies = [], dipoles_files = [], nanofq_seed = ''):
    global wdir
    global reference 
    global dip_files
    global nanofq
    wdir = workdir
    reference = reference_energies
    dip_files = dipoles_files
    nanofq = nanofq_seed
#
def fitness_function(energy,reference):
    return (1.0/np.linalg.norm(np.array(energy)-np.array(reference)))
#
#
def run_genetic_algorithm(nanofq,reference):
    #
    #
    #
    fitness_function = fq_run_and_fit
    #
    # genes = chi, eta, rq
    #
    ga_instance = pygad.GA(num_generations=20,
                           num_parents_mating=2,
                           fitness_func=fitness_function,
                           sol_per_pop=30,
                           num_genes = 4,
                           gene_space = {'low': 0, 'high': 1}
                           )
    #
    ga_instance.run()
    #
    solution = ga_instance.best_solution()[0]
    #
    new_embedding = polarizable_embedding_class.polarizable_embedding(pqeq = False)
    new_embedding.forcefield = 'fq'
    new_embedding.atomtypes = ['O','H']
    new_embedding.chi = [i for i in solution[0:2]]
    new_embedding.eta = [i for i in solution[2:4]]

    #target_directory = wdir+ 'best/'
    #os.mkdir(target_directory)
    #nanofq_interface.create_input(molecule_file = mol_dir + mol_file,   \
    #                              EEdipoles_file = dip_dir + dip_file,  \
    #                              which_dipoles = [0],                  \
    #                              nanofq_output_file_name = 'test.mfq', \
    #                              target_directory = target_directory,  \
    #                              polarizable_embedding = new_embedding,\
    #                              computation_name = 'test')
    #
    #ga_instance.plot_fitness() 
    #
    return new_embedding
#
def fq_run_and_fit(ga_instance,solution,solution_idx):
    #
    """The fitness procedure"""
    #
    # qui mi arriano le soluzioni che sono i miei geni, a questo punto, dai geni dovro' inizializzare
    # un polarizable_embedding, creare il conto e lanciarlo
    #
    new_embedding = polarizable_embedding_class.polarizable_embedding(pqeq=False)
    new_embedding.forcefield = 'fq'
    new_embedding.atomtypes = ['O','H']
    new_embedding.chi = [i for i in solution[0:2]]
    new_embedding.eta = [i for i in solution[2:4]]
    # 
    target_directory = wdir+ 'g' + str(ga_instance.generations_completed)+'_p' + str(solution_idx)
    os.mkdir(target_directory)
    energy = []
    # INIZIO CICLO
    for dip_file in dip_files:
        dipoles = dipoles_class.dipoles()
        dipoles.initialize_from_dip(dip_file)
        #
        which_dipoles = get_which_dipoles_from_dip(dip_file)
        #
        new_nanofq = nanofq_class.nanofq(molecule = nanofq.molecule, dipoles = dipoles, nanofq_path = nanofq.nanofq_path)
        #
        new_nanofq.which_dipoles = which_dipoles
        new_nanofq.polarizable_model = new_embedding
        #
        calc_name = new_nanofq.guess_name_from_dip()
        new_nanofq.name = target_directory +  '/' + calc_name
        #
        new_nanofq.create_input(input_ = new_nanofq.name + '.mfq', computation_comment = new_nanofq.name, \
                                which_dipoles = which_dipoles)
        #
        new_nanofq.run()
        #
        energy.append(new_nanofq.get_energy())
    #
    # FINE CICLO
    #
    #
    print('generation: ' + str(ga_instance.generations_completed) + ' member: ' + str(solution_idx) + \
          ' energy diff: ' + str(np.linalg.norm(np.array(energy)-np.array(reference))))
    #
    fitness = fitness_function(energy,reference)
    #
    subprocess.run(['rm', '-rf', target_directory])
    #
    return fitness
#
#
#
def get_which_dipoles_from_dip(dipfile):
    #
    #
    #
    name = dipfile.split('.dip')[0]
    name = name.split('/')[-1]
    name = name.split('dip')[1]
    name = name.split('_')[0]
    return [int(name)]
