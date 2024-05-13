![Status](https://img.shields.io/badge/Status-Work_in_Progress-yellow)
> **_NOTE:_**  This project is in progress, additional and more detailed documentation will be soon added.


## FQ parameterization
This python project is created for the Lab. sessions of the course Advanced Topics in Quantum Chemistry (Prof. T. Giovannini, A.A. 2023/2024) to help you out in the parameterziation of atomic Electronegativity (χ) and Chemical Hardness (η) for the Fluctuating Charges force field for QM/MM calculations.

## Theoretical Introduction
The [Fluctuating Charges](https://pubs.aip.org/aip/jcp/article/101/7/6141/166365/Dynamical-fluctuating-charge-force-fields) (`FQ`) model is a molecular dynamics model in which solvent atoms are endowed with a polarizable partial charge. 
The value of such partial charges depend upon atomic **Electronegativity (χ)**, **Chemical Hardness (η)**, and their inteaction with the surrounding environment (other partial charges, external potentials, ...).

The FQ model is widely employed in [Quantum Mechanics / Molecular Mechanics](https://pubs.aip.org/aip/jcp/article/157/21/214101/2842082/Assessing-the-quality-of-QM-MM-approaches-to) (`QM/MM`) models for the simulations of the spectoscopic properties of molecular systems embedded in an external environment.
Usually, the target analyte is treated at the QM level of theory employing e.g. [Density Functional Theory](https://en.wikipedia.org/wiki/Density_functional_theory), whereas the solvent environment is treated by means of the FQ model. 
By solving a set of coupled QM/MM linear equations, encoding also all the QM and MM mutual interactions, one is capable to compute the solute's excitation energies and absorption spectrum

However, there is a caveat. 
The values of **Electronegativity (χ)** and **Chemical Hardness (η)** are available only for [specific solvents](https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00763), thus limiting the applicability of FQ.

Hence the need of a user-friendly code for FQ χ, and η parameterization for any solvent.

## Methodology
To find the  **Electronegativity (χ)** and **Chemical Hardness (η)** values, we will try to reproduce some reference QM calculations (Electrostatic + Polarization contribuitions to the interaction energy, and the molecular polarizability tensor).
Such calculations are generated within the code by interfacing with [eT program](https://etprogram.org/).
Then we deploy a genetic algorithm approach to find the correct (χ,η) values matching the reference QM calculations.

## Documentation
List of the presentations of the Lab. 
Sessions of the previous year:
1. [**Session 1**](https://github.com/nicoli-luca/fq_parameterization/blob/main/docs/Session_1.pdf): Definition of the QM/FQ equations, introduction to genetic algorithms and FQ parameterization, positioning of fixed charges around target molecule;
3. [**Session 2**](https://github.com/nicoli-luca/fq_parameterization/blob/main/docs/Session_2.pdf): Definition of electrostatic embeddings (EE), QM/EE and FQ/EE working equations, setup and calculation on `HPC` of reference QM/EE calculations;
4. [**Session 3**](https://github.com/nicoli-luca/fq_parameterization/blob/main/docs/Session_3.pdf): Genetic Algorithm (GA) in the specific case of the optimization of FQ parameters. Tests using GA to compute FQ/EE calculations;
6. [**Session 4**](https://github.com/nicoli-luca/fq_parameterization/blob/main/docs/Session_4.pdf): Finalization of the GA for the optimization of FQ parameters. Parameterization of FQ for liquid water, parameterization of FQ for Deep Eutectic Solvents.

## Author
Luca Nicoli
