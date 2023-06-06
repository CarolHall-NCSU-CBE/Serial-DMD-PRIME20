# PRIME20 - Coarsed-grained force field with discontinuous molecular dynamics simulations for peptide self-assembly modelling 
## Table of contents
* [Introduction](#introduction)
* [Requirement and Installation](#requirement-and-installation)
* [Getting Started](#getting-started)
* [Running Simulation](#running-simulation)
* [Developing Status](#developing-status)
## Introduction
PRIME20 is a coarse-grained, implicit-solvent, intermediate-resolution protein model that was developed by the Hall group at North Carolina State University. The model was designed to be used with discontinuous molecular dynamics simulations (DMD) to investigate self-assembly of short peptides from their random denatured states. PRIME20 contains geometric and energetic parameters that describe the sidechain-sidechain interactions of all 20 natural amino acids. In PRIME20, each amino acid is represented by four beads: one for the amino group (NH), one for the alpha carbon (CαH), one for the carbonyl group (CO), and one for the side chain (R). DMD/RPIME20 simulation systems are canonical ensemble (NVT) with constant number of molecules (N), simulation box volume (V) and simulation temperature (T). Temperature is maintaned by using Anderson thermostat. Neutral pH water solvent is described implicitly within the force-field. Peptides that are built and simulated by PRIME20 are capped at both terminus. DMD/PRIME20 has been used successfully to simulate spontaneous α-helix, β-sheet, and amyloid fibril formation starting from the denatured conformations of peptides such as prion proteins fragments, tau protein fragments, Aβ16-22 peptides, and  Aβ17-42 peptides.

## Requirement and Installation
- The package has been developed since 2001 using Fortran90
- Fortran Intel compiler `ifort` is required.
- The installation is through the terminal.
- The source codes are in `/src/`. To compile, go to `/src/` directoy on your local device and hit `make` 

## Getting Started
**/example/**: this directory contains an example of required file and subdirectories for a simulation using DMD/PRIME20.
Requirements to start a simulation including:
- **input.txt**: Please follow the format to enter all parameters that are required for a simulation. The explanation for each parameters are also included in the file.

>Note 1: The current version only allows annealing simulation with a fixed set of temperatures. Please do not change the value of 'annealing'. Upcoming version will allow user to define annealing temperatures and time to run annealing simulation.

>Note 2: If an error is returned and the simulation is terminated during the generating of initital configuration. Adding another parameter to the end of **input.txt**: 
>
>	*sidechainmove* = value that is larger than 3.0	
>	
>It is recommended to increase only 0.5 at a time starting from 3.0. A very large number will make the initial configuration generation very slow`
- 5 empty directories for data recording must be created before submitting a job. The names of these directories must be exact.
	- `/checks/`: files for checking if the initial configuration is created correctly
	- `/inputs/`: files to record residue id and positions for each peptide sequence  
	- `/outputs/`: output files for each simulation round
	- `/parameters/`: sidechain parameters generated from the inital configuration step that are required for simulation steps
	- `/results/`:  simulation results for data analysis
		1. .bptnr: collision, bond partner of each particle
		2. .config: collision, time, particle coordinates
		3. .energy: collision, time, kinetic energy, total energy, etc.
		4. .lastvel: collision, velocities 
		5. .pdb: pdb file
		6. .rca: distance from sidechain to each particle in the backbone of a residue
>Note: These subdirectories in the **/example/** directory contains results from a short simulation for your reference. When running a new simulation, these subdirectories must be empty to avoid incorrectly data appending. When running a continueing simulation, keep all results from previous simulation in these directories. 
## Running simulation
DMD simulation using PRIME20 starts with building initial configuration. The current version is effective for system of no more than 31-residue peptides. It is recommended that concentration and number of peptide chains are reduced for longer peptides to avoid overlap due to overcrowding. User should check output file for overlapping error and reduce system size (number of peptides or concentration) if error is reported. PRIME20 allows simulations of a homogenous system or a heterogeneous system of two different peptides.

### Submit a job:
Steps to submit a simulation is as follow. These steps are after the package is succesfully installed on your device and *the path to executable file is obtained*.
1. Make a directory to run simulation or copy over the /'example/' directory, rename and then delete all files within subdirectories. If making new directory, follow the next steps. 
2. In this directory, make an 'input.txt' file following the example. You can copy over this file and change the parameters correspoding to your system.
3. In this directory, make 5 empty subdirectories at listed above if running a new simulation, or copy over these subdirectories with all data in them for a continuing simulation. 
4. Submit job. It is not recommended to run DMD/PRIME20 on terminal as a job can take days to finish. A simple bash script (.sh) to submit job is attached in '/example/'. The format is as follow.  
> #!/bin/bash
> 
> /**path_to_executive_file_DMDPRIME20**/DMDPRIME20

The bold line will need to be changed to the path to your executable file 'DMDPRIME20'. For example: If you save the package to '/home/user/Serial-DMD-PRIME20' then the path to executable file will be '/home/user/Serial-DMD-PRIME20/src/'. Your submission script will be:
> #!/bin/bash
> 
> **/home/user/Serial-DMD-PRIME20/src**/DMDPRIME20

At the beginning of DMD simulation, the system will be heated to a high temperature and then be slowly annealed to the desired temperature. This step is to make sure that all peptide chains are denatured and that the DMD simulation starts with all random coils. The numbers of collisions are defined by users. Larger system will need longer simulation times. It is recommended to start the simulation with no longer than 100 billion collisions. If the system has not aggregated after 100 billion collision, the simulations can be extended.

## Developing Status
The software is being developed and updated. An result analysis package is being developed.
