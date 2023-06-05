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
- The source codes are mainly in `/src/`. To compile, got to `/src/` directoy on your local device and hit `make` 

## Getting Started
**/example/**: this directory contains an example for a simulation using DMD/PRIME20.
Requirements to start a simulation including:
- **input.txt**: Please follow the format to enter all parameters that are required for a simulation. 
`Note: If an error is returned and the simulation is terminated during the generating of initital configuration. Adding another parameter to the end of **input.txt**: 
	sidechainmove = *value that is larger than 3.0*
It is recommended to increase only 0.5 at a time starting from 3.0. A very large number will make the initial configuration generation very slow`
- 5 directories for data recording must be created before submitting a job. The names of these directories must be exact.
	- `/checks/`: files for checking if the initial configuration is created correctly
	- `/inputs/`: files to record residue id and positions for each peptide sequence  
	- `/outpus/`: output files for each simulation round
	- `/parameters/`: sidechain parameters generated from the inital configuration step that are required for simulation steps
	- `/results/`:  simulation results for data analysis
		1. *.bptnr: collision, bond partner of each particle
		2. *.config: collision, time, particle coordinates
		3. *.energy: collision, time, kinetic energy, total energy, etc.
		4. *.lastvel: collision, velocities 
		5. *.pdb: pdb file
		6. *.rca: distance from sidechain to each particle in the backbone of a residue

## Running simulation
DMD simulation using PRIME20 starts with building initial configuration. The current version is effective for system of no more than 31-residue peptides. It is recommended that concentration and number of peptide chains are reduced for longer peptides to avoid overlap due to overcrowded. User should check output file for overlapping error and reduce system size (number of peptides or concentration) if error is reported. PRIME20 allows simulations of a homogenous system or a heterogeneous system of two different peptides.

### II. Submit a job:
*The following submission steps are writen for submitting job in Linux system. The command **nohup** submit the job to run in the background and still run the simulations when the user logs out. If using different system, the user will need to modify the script and use different submission command corresponding to their system.* 

Steps to submit a simulation is as follow:
1. Delete **nohup.out** file before the simulation to prevent new simulation output from appending to old nohup.out file.
2. Make the submission script excecutable: **chmod +x submissionscript.sh**
3. Submit job: **nohup ./submissionscript.sh &**

An example of submissionscipt.sh is as follow.

![Temp Doc/images/submissionscript.png](https://github.com/CarolHall-NCSU-CBE/Serial-DMD-PRIME20/blob/5eaa761bcdac4380ae3ee64845596951d801e78b/Temp%20Doc/images/submissionscript.png)

At the beginning of DMD simulation, the system will be heated to a high temperature and then be slowly annealed to the desired temperature. This step is to make sure that all peptide chains are denatured and that the DMD simulation starts with all random coils. This first loop in the script which is shown in the below image should not be changed for any simulation. 

![Temp Doc/images/annealing.png](https://github.com/CarolHall-NCSU-CBE/Serial-DMD-PRIME20/blob/8ebe9e46a5c20129c74ce8ccb5cc311bd75873a2/Temp%20Doc/images/annealing.png)

The numbers of collisions are defined by users. Larger system will need longer simulation times. It is recommended to start the simulation with no longer than 100 billion collisions. If the system has not aggregated after 100 billion collision, the simulations can be extended. When extend simulation time, resubmit the script with all part of the script commneted out, except the last loop in the script.

![Temp Doc/images/simulationloop.png](https://github.com/CarolHall-NCSU-CBE/Serial-DMD-PRIME20/blob/0b52f15932624b4a49c927d5baba649b843e7876/Temp%20Doc/images/simulationloop.png)

## Developing Status
The software is being developed and updated.  
