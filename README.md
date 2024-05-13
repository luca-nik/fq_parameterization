## FQ parameterization
This python project is created to help you out in the parameterziation of the Fluctuating Charges force field for QM/MM calculations.

## Theoretical Introduction
The [Fluctuating Charges](https://pubs.aip.org/aip/jcp/article/101/7/6141/166365/Dynamical-fluctuating-charge-force-fields) model is a molecular dynamics model in which solvent atoms are endowed with a polarizable partial charge. 
The value of such partial charges depend upon atomic electronegativity (χ), and Chemical Hardness (η), and their inteaction with the surroinding environment (other partial charges, external potentials, ...).

The FQ model is widely employed in [Quantum Mechanics / Molecular Mechanics](https://pubs.aip.org/aip/jcp/article/157/21/214101/2842082/Assessing-the-quality-of-QM-MM-approaches-to) (`QM/MM`) models for the simulations of the spectoscopic properties of molecular systems embedded in an external environment.
Usually, the target analyte is treated at the QM level of theory employing e.g. [Density Functional Theory](https://en.wikipedia.org/wiki/Density_functional_theory), whereas the solvent environment is treated by means of the FQ model. 
By solving a set of coupled QM/MM linear equations, encoding also all the QM and MM mutual interactions, one is capable to compute the solute's excitation energies and absorption spectrum

However, there is a caveat. 
The values of **Electronegativity (χ)** and **Chemical Hardness (η)** are available only for [specific solvents](https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00763), thus limiting the applicability of FQ. 



## Description

## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.

## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew. However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a specific context like a particular programming language version or operating system or has dependencies that have to be installed manually, also add a Requirements subsection.

## Usage
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably include in the README.

## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Luca Nicoli

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
