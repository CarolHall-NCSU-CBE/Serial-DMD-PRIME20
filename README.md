# PRIME20 - Coarsed grained force field with discontinuous molecular dynamic for peptide self-assembly modelling 
## Table of contents
* [Introduction](#introduction)
* [Technologies](#technologies)
* [Setup](#setup)
* [Running Simulation](#running-simulation)
## Introduction
PRIME20 is a coarse-grained (CG) implicit-solvent intermediate-resolution protein model that was developed by the Hall group at North Carolina State University. The model was designed to be used with discontinous molecular dynamics simulations. PRIME20 contains geometric and energetic parameters that describe the sidechain-sidechain interactions of all 20 natural amino acids. In PRIME20, each amino acid is represented by four beads: one for the amino group (NH), one for the alpha carbon (CαH), one for the carbonyl group (CO), and one for the side chain (R). DMD/PRIME20 has been used successfully to simulate spontaneous α-helix, β-sheet, and amyloid fibril formation starting from the denatured conformations of peptides such as prion proteins fragments, tau protein fragments, Aβ16-22 peptides, and  Aβ17-42 peptides.
## Technologies
The package has been developed since 2001 using Fortran90
## Setup
The package doesn't require installation on your devide but the whole package must be copied over to a new directory to start a new simulation.
## Running simulation
All files that are required for DMD/PRIME20 simulation are acccesible from the directory **submissionfiles**. The files include:
1. *Inputfile.f90*
2. *Submissionscript.sh*
### I.	Generating initial configuration
DMD simulation using PRIME20 starts with building initial configuration. The current version is effective for system of less than 31-residue peptides. It is recommended that concentration and number of peptide chains are reduced for longer peptides to avoid overlap due to overcrowded. User should check output file for overlapping error and reduce system size (number of peptides or concentration) if error is reported.
The inputfile.f90 contains all the parameters that are required for a simulation. PRIME20 allows simulations of a homogenous system or a heterogeneous system of two different peptides. 
 ![Temp Doc/images/initial_allinone.png](https://github.com/CarolHall-NCSU-CBE/Serial-DMD-PRIME20/blob/45eb102c71d57b322d413f7297eed412a19df235/Temp%20Doc/images/initial_allinone.png)
1. Specified the peptides for the simulations. If the simulation system is homogeneous, parameters pep1 and pep2 are the same. 
2. Specified the number of beads within a peptide (nb1 and nb2). As PRIME20 is a 4 beads coarse-grained model, the number of beads is equal to the chain length multiplied by 4.
3. Specified the number of beads in a peptide without glycines (numbeads1 and numbead2). As glycine does not have a sidechain, numbeads1 is equal to nb1 minus the number of glycines in that peptide. Numbeads2 is found as similar.
4. Specified chain length (chnln1 and chnln2)
5. Specified how many peptide chains for each peptide (nc and nc2) 
6. Specified the length of the simulation box in Angstrom (boxlength)
boxlength=((Total number of peptide chains*1000)/(Avogadro^' s number*Concentration))^(1/3)*10^9
where: Concentration is in mM

7. Specified simulation temperature in Kelvin (simtemp)
8. The two parameters dadjust1 and dadjust2 are not recommended to be changed unless an error is returned and the simulation is terminated during the generating of initital configuration. If seeing error, slightly increase dadjust1 and dadjust2. An example of the error is as follow.
 ![Temp Doc/images/Error.png]
Results of the initial configuration generation step are recorded in different sub_directories
•	Within genconfig directory:
	-Compiled file
	-Output file
	-chninfo-n1.data
	-chninfo-n2.data
Notes: these files need to be deleted before the new system is generated if the entire package is copied over.
•	Within the sub_directory results:
	-run0000.lastvel
	-run0000.energy
	-run0000.config
	-run0000.bptnr
Notes: These files contain initial configurations and velocities
•	Within the sub_directory parameters:
	-identity.inp
	-hp1.inp
	-hp2.inp
	-firstside1.data
	-firstside2.data
Notes: These files contain identity of peptides in the system
•	the sub_directory check: All files in here are for the users to check configuration, velocities, mass and energy of the initial system.


### II. Submit a job:
To start the simulation, submit the bash script submissionscript.sh using the following command
	nonhup ./submissionscript.sh &
Note: before submit the script, make sure that the script is executable chmod +x submissionscript.sh\
An example of submissionscipt.sh is as follow. The script is written for linux system.
 
At the beginning of DMD simulation, the system will be heated to a high temperature and then be slowly annealed to the desired temperature. This step is to make sure that all peptide chains are denatured and that the DMD simulation starts with all random coils. This highlighted section in the below image should not be changed for any simulation.
 
The numbers of collisions are defined by users. Larger system will need longer simulation times. It is recommended to start the simulation with no longer than 100 billion collisions. If the system has not aggregated after 100 billion collision, the simulations can be extended. When extend simulation time, resubmit the script with all part of the script commneted out, except the last loop in the script.



### III. Results:
All results are saved in the results directory and can be read for data analysis:
	*.config: collision, time, particle coordinates
	*.lastvel: collision, velocities 
	*.bptnr: collision, bond partner of each particle
	*.energy: collision, time, kinetic energy, total energy, etc.
	*.pdb

