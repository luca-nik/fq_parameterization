from classes import molecule_class
from classes import dipoles_class
from classes import polarizable_embedding_class
import constants
import numpy as np
import sys
import os
#
def run_algorithm():
    #
    #
    # 1) Avere i conti et fatti
    # 2) Inizializzare i conti nanofq (o avere una cartella iniziale)
    # 3) Fare un ciclo in cui lanciamo i conti nano_fq
    # 4) Estrarre le energie e calcolare la loss_function
    #https://pygad.readthedocs.io/en/latest/README_pygad_ReadTheDocs.html#stop-at-any-generation

