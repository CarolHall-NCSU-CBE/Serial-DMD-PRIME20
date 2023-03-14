# PRIME20 - Coarsed-grained force field with discontinuous molecular dynamics simulations for peptide self-assembly modelling 
## Table of contents
* [Introduction](#introduction)
* [Requirement](#requirement)
* [Setup](#setup)
* [Running Simulation](#running-simulation)
* [Developing Status](#developing-status)
## Introduction
PRIME20 is a coarse-grained, implicit-solvent, intermediate-resolution protein model that was developed by the Hall group at North Carolina State University. The model was designed to be used with discontinuous molecular dynamics simulations (DMD) to investigate self-assembly of short peptides from their random denatured states. PRIME20 contains geometric and energetic parameters that describe the sidechain-sidechain interactions of all 20 natural amino acids. In PRIME20, each amino acid is represented by four beads: one for the amino group (NH), one for the alpha carbon (CαH), one for the carbonyl group (CO), and one for the side chain (R). DMD/RPIME20 simulation systems are canonical ensemble (NVT) with constant number of molecules (N), simulation box volume (V) and simulation temperature (T). Temperature is maintaned by using Anderson thermostat. Neutral pH water solvent is described implicitly within the force-field. Peptides that are built and simulated by PRIME20 are capped at both terminus. DMD/PRIME20 has been used successfully to simulate spontaneous α-helix, β-sheet, and amyloid fibril formation starting from the denatured conformations of peptides such as prion proteins fragments, tau protein fragments, Aβ16-22 peptides, and  Aβ17-42 peptides.

## Requirement
- The package has been developed since 2001 using Fortran90
- Fortran Intel compiler ifort is required.
- Submission script is writen for Linux and will require modification if running on different system.

## Setup
- For every new simulation, the whole package must be downloaded, cloned, or copied over to the local machine where the simulation is run.

## Running simulation
Notes: Currently, the user need to access the source codes to specify parameters for simulation as well as to analyse data. Newer version will be soon updated. 
- For current version, all files that are required for DMD/PRIME20 simulation are acccesible from the directory **submissionfiles**. The files include:
1. *Inputfile.f90*
2. *Submissionscript.sh*
- All results are saved in the results directory and can be read for data analysis:
1. *.config: collision, time, particle coordinates
2. *.lastvel: collision, velocities 
3. *.bptnr: collision, bond partner of each particle
4. *.energy: collision, time, kinetic energy, total energy, etc.
5. *.pdb

### I.	Generating initial configuration
DMD simulation using PRIME20 starts with building initial configuration. The current version is effective for system of no more than 31-residue peptides. It is recommended that concentration and number of peptide chains are reduced for longer peptides to avoid overlap due to overcrowded. User should check output file for overlapping error and reduce system size (number of peptides or concentration) if error is reported.
The file *inputfile.f90* contains all the parameters that are required for a simulation. PRIME20 allows simulations of a homogenous system or a heterogeneous system of two different peptides. The following image shows *inputfile.f90* file.

 ![Temp Doc/images/initial_allinone.png](https://github.com/CarolHall-NCSU-CBE/Serial-DMD-PRIME20/blob/45eb102c71d57b322d413f7297eed412a19df235/Temp%20Doc/images/initial_allinone.png)
 
Users must specify all required input parameters below corresponding to the system they want to simulate.
1. **pep1**: peptide sequence for component 1 in one letter abbreviating format (*eg. pep1 = "GVLYVGS"*)
2. **pep2**: peptide sequence for component 2 in one letter abbreviating format. If the system is homogeneous, pep2 is same as pep1.
3. **nb1** and **nb2**: number of beads of peptide 1 and 2, respectively. As PRIME20 is a 4 beads coarse-grained model, the number of beads is equal to the number of residues multiplied by 4.
4. **numbeads1** and **numbead2**: number of beads in peptide 1 and peptide 2 without glycines, respectively. As glycine does not have a sidechain, 

	$$ **numbeads1** = {**nb1** - {number\ of\ glycines\ in\ that\ peptide}} $$
	
	$$ **numbeads2** = {**nb2** - {number\ of\ glycines\ in\ that\ peptide}} $$	

4. **chnln1** and **chnln2**: number of residues in peptide 1 and peptide 2, respectively.
5. **nc** and **nc2**: number of peptide chains for each peptide. If the system is homogeneous, it's recommended that nc = nc2 = total peptide chains divides by 2.
6. **boxlength**: length of the simulation box in Angstrom


```
code
```
$$ boxlength = \left\lbrack{\\frac{{Total\ number\ of\ peptide\ chains} * {1000}}{{Avogadro's\ number} * {Concentration}} }\right\rbrack^{1/3} * 10^9 $$

	where: Concentration is in mM

7. **simtemp**: simulation temperature in Kelvin
8. The two parameters **dadjust1** and **dadjust2** are not recommended to be changed unless an error is returned and the simulation is terminated during the generating of initital configuration. If seeing error, slightly increase **dadjust1** and **dadjust2** by *0.5 at a time*. An example of the error is as follow.

 ![Temp Doc/images/Error.png](https://github.com/CarolHall-NCSU-CBE/Serial-DMD-PRIME20/blob/ace39b9324962999c9f1ee448907000c8d65d9e1/Temp%20Doc/images/Error.png)
 
DMD/PRIME20 builts simulation system by generating the first peptide of each type and then replicates those first peptides randomly in the periodic boundary condition box with the specified length.
 
Results of the initial configuration generation step are recorded in different sub_directories
- Within genconfig directory:
	- Excecutable file: the configuration generation step
	- chninfo-n1.data: coordinates of all beads in the first peptide 1  
	- chninfo-n2.data: coordinates of all beads in the first peptide 2
- Within the sub_directory results:
	- run0000.lastvel: initial velocity of all beads in the system
	- run0000.config: initial coordinates of all beads in the system
	- run0000.energy: energy file that is empty before the simulation starts
	- run0000.bptnr: bond partner of each bead in the system that is empty before the simulation starts
- Within the sub_directory parameters:
	- identity.inp: bead IDs. In DMD/PRIME20, each bead type (NH, CO, Calpha, or 20 sidechains) is assigned a numerical ID. The bead names, bead IDs and bead masses are documented in *parameters/mass.data*.
	- hp1.inp: hydrophobicity of each beads in peptide 1. 1 means hydrophobic and 0 means non-hydrophobic.
	- hp2.inp: hydrophobicity of each beads in peptide 2. 1 means hydrophobic and 0 means non-hydrophobic.
	- firstside1.data: locating glycine in the peptide 1. 0 means glycine, 1 means all other peptides. 
	- firstside2.data: locating glycine in the peptide 2. 0 means glycine, 1 means all other peptides.
Notes: These files contain identity of peptides in the system
- the sub_directory check: All files in here are for the users to check configuration, velocities, mass and energy of the initial system.

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
