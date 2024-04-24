# Serial DMD/PRIME20 - Coarse-grained force field for use with discontinuous molecular dynamics simulations to model peptide self-assembly 
## Table of contents
* [Introduction](#introduction)
* [Methods](#methods)

	[Coarse-grained Model](#coarse-grained-model)

 	[DMD Simulation](#dmd-simulation)
* [Requirement and Installation](#requirement-and-installation)
* [Units](#units)
* [Running Simulation](#running-simulation)

  	[Getting Started](#getting-started)
  
  	[File and folder requirements](#file-and-folder-requirements)
  
  	[Paramemters in input.txtx to start a simulation](#paramemters-in-input.txtx-to-run-a-simulation)

  	[Submit A Job](#submit-a-job)
* [Analysis Package](#Analysis-Package)
* [Developing Status](#developing-status)
## Introduction
PRIME20 is a implicit-solvent and intermediate-resolution protein model that was developed by the Hall group at North Carolina State University. The model is designed to use with discontinuous molecular dynamics (DMD) simulations to investigate self-assembly of short peptides from their random denatured states. PRIME20 contains pair-wise interaction parameters of all twenty natural amino acid. In PRIME20, each amino acid is represented by four beads: one for the amino group (NH), one for the alpha carbon (CαH), one for the carbonyl group (CO), and one for the side chain (R). DMD/PRIME20 simulation systems are canonical ensemble (NVT) with constant number of molecules (N), simulation box volume (V) and simulation temperature (T). Temperature is maintaned by using Anderson thermostat. Neutral pH water solvent is described implicitly within the force-field. Peptides that are built and simulated by PRIME20 are capped at both terminus. DMD/PRIME20 has been used successfully to simulate spontaneous α-helix, β-sheet, and amyloid fibril formation starting from the denatured conformations of peptides such as prion proteins fragments[1][2], tau protein fragments[3], Aβ16-22 peptides[4][5][6], and  Aβ17-36 peptides[7]. DMD/PRIME20 is also used with PepAD - a peptide design software, to discover amyloid-forming peptides [8].

## Methods
### Coarse-grained model:
PRIME20 is a four-bead CG model with three-bead to describe backbone and one-bead to describe sidechain. The three-bead backbone model is found to be a well balance between relatively realistic backbone geometry and model efficiency. Firstly, it can explicitly model hydrogen bonding in the backbone. Secondly, it allows modeling of secondary structure of polypeptide as dihedral angles can be described, calculated, and constrained in the sterically allowed spaced on the Ramachandran plot. The backbone beads are assigned different diameters and mass based on the realistic size of the groups they are presented in the model. The amine bead is positioned at the corresponding nitrogen atom, while carbon alpha and carbonyl bead are positioned at the respective carbon atom of that group. While the backbone is modeled with details, each sidechain is mapped into a single united bead that is positioned at the center of mass of the sidechain group. Sidechain masses are kept as realistic values, but the diameter of each sidechain is parameterized specifically for each sidechain-sidechain interaction. Cheon et al.[9] calculated radial distribution function histogram between each pair of sidechains from 711 pdb of natural proteins from Protein Data Bank and selected the first non-zero value of the radial distribution function to be the diameter of the pair. Cheon et al. tabulated a set of 171 bead diameters as glycine doesn’t have a sidechain. In summary, the diameter of the sidechain is not a constant value through the simulation but will be changed based on the interaction sidechain. Covalent and pseudo bonds are implemented in the models to secure modeled polypeptides in realistic geometry of natural proteins. Three covalent bonds connect NH, CO, and R to CαH, and one covalent bond describes peptide bond linking two peptides together. These bonds are shown as solid black bonds in Fig. 1. The lengths of these bonds follow the corresponding distances in real proteins. Instead of using a set of bond angles, six pseudo bonds are assigned for both inter- and intra-residues to manage the peptide conformations as L-isomer and dihedral bond angles in the realistic observed regions. The pseudo bond that connects CαH-CαH is set at 3.8 Angstroms to restrain the trans configuration of the polypeptide. Currently, PRIME20 doesn’t consider cis configurations, therefore PRIME20 cannot model proline as accurately as other amino acids. Numeric values of these bond lengths are listed by Voegler, Smith and Hall.
### DMD Simulation:
DMD is very fast compared to classical MD as DMD uses simple discontinuous potentials such as hard sphere or square-well potentials to describe the interaction between CG beads. DMD is event-driven as opposed to molecular dynamics which is time-driven. An event or also called a collision is said to occur whenever there is an abrupt change in the interaction potential between two beads. In this paper, we use the term collision as one of the time units beside reduced time and time SI unit. DMD simulations of polymers and other chain-like molecules were first introduced by Rapaport[10-12] and later modified by Bellemans et al.[13] DMD simulations are later used by man research groups to study protein aggregation. Each group uses DMD with different forcefields, including CG and all-atom.
DMD speeds up simulation by both increasing timestep and reducing the cost of calculations. At each iteration, the algorithm assumes that all beads in the system move at their constant velocities meaning that there is no force applied on the beads during movement. Therefore, the DMD method reduces computational cost to numerically solve the Newton’s equation of motion at every small timestep (1-2fs). Instead, the program solves for collision time of each pair when the distance between them is at the discontinuity of the discontinuous interaction potential. Single-event scheduling is used to store the collision time. Therefore, the algorithm only stores the soonest-to-happen event for each bead. Bucket sorting algorithm is used to reduce the cost of searching for the soonest-to-happen event for the whole system, as only a small portion of collision times that are stored in the first bucket goes through bubble sorting. DMD then solves dynamic changes of the two collided beads using conservation of energy and momentum. The system is advanced to that soonest-to-happen collision time.
Different types of events might occur including excluded-volume events, bond events, or square-well events based on types of collided beads and direction of their velocities. Excluded-volume events are observed when an infinite repulsion between adjacent beads occurs if they move too close to each other. Bond events occur when bonded beads move outward and try to break their assigned bonds. As a result, infinitely strong attractive forces are applied so that the beads move back and forth within permissible bond lengths. Square-well events include more than just one type of interaction between two beads. First is the well-capture event when a sphere enters the square well of another sphere. Another is well-bounce, when a sphere attempts to leave the potential well but, due to insufficient kinetic energy, it is unable to overcome the negative potential energy of the attractive square well and is reflected into the well. Well-dissociation occurs when the bead has sufficient kinetic energy to completely leave the square well of the other sphere, resulting in a decreased kinetic energy of the bead. 

## Requirement and Installation
- Requirement: Fortran F90 and Intel Fortran Compiler ifort
- The installation is through the terminal.
- The source codes are in /src/. To compile, open a terminal and then navigate to the /src/ directory on your local device. Once in /src/ directory, create the executed files by enter the commands below.
- To create `initconfig` for generating initial configuration
>
	make -f genconfig.mk

- To create `DMDPRIME20` for DMD simulations
>
	make -f dmd.mk

- To create `DMDanalysis` for data anylysis
>
	make -f dmd_analysis.mk

If there is no error return, check if `initconfig`, `DMDPRIME20`, and `DMDanalysis` are succesfully created in src
Obtain the paths to these executable files to use in job submission.
>**Note:** if redownload the package or update a new version, the previous steps need to be redo.

## Units
Table 1: Units that are used in `input.txt` and result analysis  
|Quantity   |Unit                                                      |
|-----------|----------------------------------------------------------|
|boxlength  | Angstrom                                                 |
|collision  | collision in `input.txt` or billion collisions in results|
|energy     | kJ/mol                                                   |
|temperature| K                                                        |
|time       | microsecond                                              |


## Running simulation
### Getting Started:
#### File and folder requirements
- All files and folders that are required for a complete simulation and data analysis can be found in directory **/example/** in the package. Names of files and folders need cannot be modified. The **input.txt**, **submission_script**. and required folders must be placed in the same directory for each individual simulation. 
- Requirements for running simulation: **input.txt** and **submission_script**.
> Format of **input.txt** must be followed exactly. **submission_script** should be written to suit the local device.
- Requirements for data recording: 5 empty directories to record simulation ouputs - `/checks/`, `/inputs/`, `/outputs/`, `/parameters/`, and `/results/`. These folders must be created at the beginning of a new simulation.
- Requirement for data analysis: 1 empty directory to record data collected from using data analysis package - `analysis`

#### Paramemters in `input.txtx` to start a simulation
**input.txt**: Please follow the format to enter all parameters that are required for a simulation. The explanation for each parameters are also included in the file.
Table 2: Paremeters for DMD/PRIME20 simulation
|Parameter            | Description                                                              |
|**pep1** and **pep2**| sequences of the peptides that are simulated. It must be in abbrevating alphabetical format (e.g. pep1=GVLYVGS).| 
|		      |The current version can run simulations for system with single or double components; each with maximum length of |
|		      |30 residues. If system contains single peptide sequence, then *pep1* and *pep2* are the same in the 'input.txt'  |

	- **chain1** and **chain2** are the number of peptide chains of each peptide component in the simulation box. If the peptide is long, *chain1* and *chain2* should be reduced to avoid overcrowding, overlapping and to reduce simulation time. The largest system has been simulated using DMD/PRIME20 contains 200 peptides chains.

	- **boxlength** is the length of the simulation box. DMD/PRIME20 uses cubic box with periodic boundary condition for all simulations. *boxlength* is selected based on the number of peptide chains and concentration:

$$ boxlength = (\frac{\text{Total number of peptide chains}*1000}{\text{Avogadro's number * Concentration}})^\frac{1}{3}*10^9 $$
- where *Concentration* is in *mM* and *boxlength* is in *Angstrom*

	- **T** is simulation temperature in *Kelvin*. When start simulations for a new system, it is recommended to run multiple simulations of the same system at different temperatures. Check the simulation results to select the temperature that predict high order peptide aggregation. The simulation might get stuck in local miminima if the temperature is too low, but there is no aggregation if the temperature is too low.
 	- **coll** which is the number of collisions for DMD/PRIME20 to finish a *round* and record simulation results. DMD/PRIME20 is designed to run, complete and record in many rounds to avoid large result files and to allow the simulation to restart if it is crashed midway. As DMD is discontinous molecular dynamics simulation, collsion (coll) is used instead of timestep. Collision will be converted to real time when running data analysis package (underdevelopment and will be updated soon). There is not a fix value in real time for a collision.
  	- **trajrecord** which is the number of collisions for DMD/PRIME20 to record a frame of trajectory 	

	- **Annealing**: The current version allows annealing simulation with a default set of temperatures (annealing = 0) or a user-defined temperatures (annealing = 1). If using user-defined temperature, include addtional parameters below the annealing line:
 		- **startingtemp**: starting temperature for the annealing process (in *Kelvin*)
		- **endingtemp**: ending temperature for the annealing process (in *Kelvin*)
		- **tempstep**: temperature drop for each annealing cycle (in *Kelvin*)
  		- **annealingcoll**: number of collisions for each temperature in annealing process. Recommended value is from 100 million to 250 million collisions   

>**Note:** If an error is returned and the simulation is terminated during the generating of initital configuration. Adding another parameter to the very end of **input.txt**: 
>
>	*sidechainmove* = value that is larger than 3.0	
>	
>It is recommended to increase only 0.5 at a time starting from 3.0. A very large number will make the initial configuration generation very slow`

- An example of **input.txt** that include parameters for are use-defined annealing temperature is below. If running simulaiton with default annealing temperature, set annealing = 0 and delete all parameter below that line.
>
	Peptide sequence 1

	pep1=GVLYVGS
 
	Number of peptide 1 chains in the system
 
	chain1=3
 
	Peptide sequence 2
 
	pep2=GVLYVGS
 
	Number of peptide 2 chains in the system
	
 	chains=3
	
 	Box length in Angstrom

 	boxlength=159.0D0
	
 	Simulation temperature in Kelvin
	
 	T=310.0D0
	
 	Result recording frequency in collisions
	
 	coll=1000000000

   	Trajectory recording frequency

    	trajrecord = 100000
	
 	Annealing process only for new simulation: 0 is the default, 1 is the specified temperatures in Kelvin
	
 	annealing = 1
	
 	startingtemp = 1000.0
	
 	endingtemp = 375.0
	
 	tempstep = 125
	
 	annealingcoll = 100000000
  
2. **submitscript.csh** is an example of the tcsh script that is used to submit a job on an HPC system. This file will need to be modified according to users' computer system. Main content of the script is the three steps of a simulation.
>
	#Generate initial configuration for new simulation:

	/path_to_initconfig/initconfig

	#Annealing:

	foreach i (`seq 1 annealingrounds`)

		do /path_to_DMDPRIME20/DMDPRIME20 < inputs/annealtemp_$i > outputs/out_annealtemp_$i

	end

	#DMD simulations

	foreach i (`seq start end`)

		do /path_to_DMDPRIME20/DMDPRIME20 < inputs/simtemp > outputs/out_simtemp_$i

	end

- Generating initial configuration for new simulation: This step is to create a cubic box that contents the number of peptide chains defined by users, position and velocity of each particles. Outputs of this step are saved in `/inputs/`, `/parameters/`, and `/results/` directories. In `/results/`, output files from generating inital configuration are named with *0000*. These files are required for any DMD/PRIME20 simulation and need to be available in their designated locations. If restarting or resuming a simulation, this step is skipped as long as the initial configuration files are available. Although, DMD/PRIME20 simulation is parallelized, this step is done in *serial*. The path to the executable file `initconfig` must be specifed. For example: If you save the package to `/home/user/Parallel-DMD-PRIME20` then the path to executable file will be `/home/user/Parallel-DMD-PRIME20/src/`. Your submission script will look like: `/home/user/Parallel-DMD-PRIME20/src**/initconfig`
- Annealing: This step is to heat up the initial system to very high temperature and then slowly cool it down to near simulation temperature. This step is only required for simulation of a completely new system. The purpose of this step is to make sure all peptide chains are denatured and simulation starts with all random coils. There are two options for annealing:
	- Default annealing (annealing = 0 in input.txt): The annealing process will be done with a default set of temperatures. These temperatures are used in many simulations since the software was developed. If using default annealing, set **annealingrounds** in the parallelscript.sch to **9**. This means the anneanling process runs at 9 different temperatures. The temperatures and number of collision at each temperature can be found in */inputs/* directory.
 	- User-defined annealing (annealing = 1 in input.txt): The annealing process will be done with the temperature range and number of collision that are defined by user. If using this option, the **annealingrounds** is found as:

$$ \text{annealingrounds} = \frac{\text{startingtemp - endingtemp}}{\text{tempstep}} + 1 $$

- DMD Simulation: This is the DMD simulation step. The number of simulation rounds must be specified. Each simulation round will be run for a number of **coll** specied in *input.txt*. 
	- Starting a new simulation: **start** = 1 and **end** = number of simulation rounds. The total simulation time is equal **coll** times **end**, so the value of **end** is dependent on how long user wants to run the simulation for and how often user wants to record outputs. It is recomended to run simulations for about 100 billion collisions first and then extending simulation times if aggregation has not happened. For example, if running for 100 simulation rounds then set ``foreach i (`seq 1 100`)``
 	- Continuing simulation or restarting crashed simulation: **start** = (the last completed simulation round + 1) and **end** = number of simulation rounds.
  		- For countinuing simulation: if the previous simulation ends after 100 simulation rounds, then to countinue the simulation to 200 simulation collisions set ``foreach i (`seq 101 200`)``
		- For crashed or incomplete simulation: if the previous simulation was set for 100 simulation rounds, but for some reasons the simulation partilly finishes at 80 simulation rounds. User must delete all the output files relating to the incomplete simulation in */results/* and */outputs/* meaning that the last complete simulation is at the round 79. Then simultion can be restart by setting ``foreach i (`seq 80 100`)``. When restarting a simulation, generating initial configuration and annealing must be skipped.
- Both Annealing and DMD simulation are designed to utilize the benefit of parallel performance. Therefore, both commands are executed using `mpirun`. Both steps are computed using the executable file `DMDPRIME20`, so the path to this file is required to be specified similar to the example for `initconfig` Temperatures and number of collisions at each temperatures must be access from `/inputs/` directory and the output files will be saved in `/outputs/` directory. The names and locations of these files are designated and cannot be changed. 
 
3. 5 empty directories for data recording must be created before submitting a job. The names of these directories must be exact.
	- `/checks/`: files for checking if the initial configuration is created correctly
	- `/inputs/`: files to record residue id (identity.inp and identity2.inp), positions for each peptide sequence (chninfo-n1.data and chninfo-n2.data), reduced annealing temperatures (annealtemp_*), and reduced simulation temperature (simtemp)   
	- `/outputs/`: output files for each simulation round
	- `/parameters/`: sidechain parameters generated from the inital configuration step that are required for simulation steps
	- `/results/`:  simulation results for data analysis
		a. .bptnr: collision, bond partner of each particle
		b. .config: collision, time, particle coordinates
		c. .energy: collision, time, kinetic energy, total energy, etc.
		d. .lastvel: collision, velocities 
		e. .pdb: pdb file
		g. .xyz: trajectory files
		f. .rca: distance from sidechain to each particle in the backbone of a residue

>Note: These subdirectories in the **/example/** directory contains results from a short simulation for your reference. When running a new simulation, these subdirectories must be empty to avoid incorrectly data appending. When running a continuing simulation, keep all results from previous simulation in these directories. The **.out** file shows an example of successful initial configuration generation. If your screen-written output look like this and no error showed, the initial configuration is successulffy generated. This *.out* file must be deleted before any simulation if it exists to avoid being confused by old data.

### Submit a job:
Steps to submit a simulation is as follow. These steps are after the package is succesfully installed on your device and *the path to executable file is obtained*.
1. Make a directory to run simulation or copy over the `/example/` directory, rename and then delete all files within subdirectories and `*.out`. If making new directory, follow the next steps. 
2. In this directory, make an 'input.txt' file following the example. You can copy over this file and change the parameters correspoding to your system.
3. In this directory, make 5 empty subdirectories at listed above if running a new simulation, or copy over these subdirectories with all data in them for a continuing simulation. 
4. Submit job. It is not recommended to run DMD/PRIME20 on a login node as a job can take days to finish. A simple tcsh script (.csh) to submit job is attached in `/example/`.

## Analysis Package
DMD-PRIME20 allows running simulations for hundred of microseconds. A simulation can take a month to complete and generates a big set of data as results. Therefore, most of data is written in binary and require extra step to extract data for specific analysis. The analysis package is included implemented to allow user to allow user access data of their choices. It's currently developed, new functions will be added based on user feedback.



After executable `DMDAnalysis` is created, the package can be used from running directory in the same location at input.txt file. Although analysis package can be run on the interactive termital, it is recommended to submit analysis job in the background or using a queuing system. It's can take from few seconds to few minutes to complete the task depending on how big result data set is. For a long simulation that take approximatly a month, result data set can be very big and required extra time for the analysis package to extract requested data.

### Using analysis package:
Directory `analysis` must be created in the running directory (same location as the `results` and `output` directories) to store result files. Except the trajectory file that is in `.xyz` format (`.xtc` in future), all other data will be imported in `.txt` files and can be read by any tools of your choice. 
### Energy vs time:
>
	path_to_DMDanalysis/DMDanalysis evst start_file end_file

A file that is named in the format of `evst//start_file//to//end_file//.txt` is created in the `analysis` directory. The file has 6 collumns in the order of: collisions, time, temperature, total energy (Etotal), kinectic energy (KE) and interaction potential enerrgy (PE). As PRIME20 simulation process includes 2 steps - anneanling and simulation, users have choice to include or exclude the annealing process in the report. In the bash script used for job submition, users have defined the number of annealing rounds so they can choose the value for start_file the analysis from `0001` to include annealing or from `annealinground+1` to exclude the annealing process. The end_file can be anyfile that user wish to stop the analysis at. 


### Hydrogen bonding vs time:
>
	path_to_DMDanalysis/DMDanalysis hbvst start_file end_file

A file that is named in the format of `hbvst//start_file//to//end_file//.txt` is created in the `analysis` directory. The file has 4 columns in the order of: collisions, time, total hydrogen bonds and interpeptide hydrogen bonds formed during the simulations. Choice of start_file and end_file is similar to energy vs time.

### Pdb file at any recorded time:
>
	path_to_DMDanalysis/DMDanalysis pdb number_of_collisions_in_billions

DMD/PRIME20 allows simulation of hundered of microseconds, therfore, the results would be very large if we recorded data at every steps. In addition, bond vibrations happen more often than the main events, therefore, there is not much to observed if recorded too often. PRIME20 records coordinate of the whole system corresponding to the recording of energy and hydrogen bonding. Therefore, it is required that the input of `number_of_collisions_in_billions` must be exactly the same as the collision that is recorded in either `evst` or `hbvst` 

## Developing Status
The software is being developed and updated. An result analysis package will be updated soon.

## References:
[1] Wang, Y., Shao, Q., and Hall, C. K. N-terminal Prion Protein Peptides (PrP(120–144)) Form Parallel In-register β-Sheets via Multiple Nucleation-dependent Pathways. Journal of Biological Chemistry. Vol. 292, Issue 50. (2016). https://doi.org/10.1074/jbc.M116.744573

[2] Wang, Y., Shao, Q., and Hall, C. K. Seeding and cross-seeding fibrillation of N-terminal prion protein peptides PrP(120–144). Protein Science (2018). https://doi.org/10.1002/pro.3421

[3] Cheon, M., Chang, I., and Hall, C. K. Influence of temperature on formation of perfect tau fragment fibrils using PRIME20/DMD simulations. Protein Science (2012). https://doi.org/10.1002/pro.2141

[4] Cheon, M., Chang, I., and Hall, C. K. Spontaneous Formation of Twisted Ab16-22 Fibrils in Large-Scale Molecular-Dynamics Simulations. Biophysical Journal. Vol. 101, 2493-2501 (2011).  https://doi.org/10.1016%2Fj.bpj.2011.08.042

[5] Samuel J. Bunce et al., Molecular insights into the surface-catalyzed secondary nucleation of amyloid-β40 (Aβ40) by the peptide fragment Aβ16–22.Sci. Adv.5,eaav8216(2019).DOI:10.1126/sciadv.aav8216

[6] Yiming Wang et al., Thermodynamic phase diagram of amyloid-β (16–22) peptide. Proceedings of the National Academy of Sciences. Vol. 116, Issue 6, 2091-2096 (2019). doi:10.1073/pnas.1819592116

[7] Wang, Y., Latshaw, D. C., and Hall, C. K. Aggregation of Aβ(17–36) in the Presence of Naturally Occurring Phenolic Inhibitors Using Coarse-Grained Simulations. Journal of Molecular Biology. Volume 429, Issue 24, 3893-3908 (2017). https://doi.org/10.1016/j.jmb.2017.10.006.

[8] Xingqing Xiao and others, Sequence patterns and signatures: Computational and experimental discovery of amyloid-forming peptides, PNAS Nexus, Volume 1, Issue 5, November 2022, pgac263, https://doi.org/10.1093/pnasnexus/pgac263

[9] Cheon, M., Chang, I., Hall, C. K. Extending the PRIME Model for Protein Aggregation to All 20 Amino Acids. Proteins Struct. Funct. Bioinforma. Vol. 78, Issue 14, 2950–2960 (2010). https://doi.org/10.1002/prot.22817.

[10] Rapaport, D. C. Molecular Dynamics Simulation of Polymer Chains with Excluded Volume. J. Phys. Math. Gen. 1978, 11 (8), L213–L217. https://doi.org/10.1088/0305-4470/11/8/008.

[11] Rapaport, D. C. Molecular Dynamics Study of a Polymer Chain in Solution. J. Chem. Phys. 1979, 71 (8), 3299–3303. https://doi.org/10.1063/1.438770.

[12] Rapaport, D. C. The Event Scheduling Problem in Molecular Dynamic Simulation. J. Comput. Phys. 1980, 34 (2), 184–201. https://doi.org/10.1016/0021-9991(80)90104-7.

[13] Bellemans, A.; Orban, J.; Van Belle, D. Molecular Dynamics of Rigid and Non-Rigid Necklaces of Hard Discs. Mol. Phys. 1980, 39 (3), 781–782. https://doi.org/10.1080/00268978000100671.
