> **_NOTE:_**  This project is on progress, additional and more detailed documentation will be soon added.


## FQ parameterization
This python project is created for the Lab. sessions of the course Advanced Topics in Quantum Chemistry (Prof. T. Giovannini, A.A. 2023/2024) to help you out in the parameterziation of atomic **Electronegativity (χ)** and **Chemical Hardness (η)** for the Fluctuating Charges force field for QM/MM calculations.

## Theoretical Introduction
The [Fluctuating Charges](https://pubs.aip.org/aip/jcp/article/101/7/6141/166365/Dynamical-fluctuating-charge-force-fields) model is a molecular dynamics model in which solvent atoms are endowed with a polarizable partial charge. 
The value of such partial charges depend upon atomic **Electronegativity (χ)**, **Chemical Hardness (η)**, and their inteaction with the surroinding environment (other partial charges, external potentials, ...).

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
List of the presentations of the Lab. Sessions of the previous year:
1. [**Session 1**](.docs/Session_1) 
3. [**Session 2**](.docs/Session_2) 
4. [**Session 3**](.docs/Session_3) 
5. [**Session 4**](.docs/Session_4) 

## Authors and acknowledgment
Luca Nicoli

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
