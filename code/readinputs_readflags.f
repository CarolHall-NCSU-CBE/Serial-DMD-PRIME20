#include "readflags.f90"
subroutine readinputs
	use inputreadin
	use global
	implicit none
	real redsimT, annealtemp(9)
	integer i, annealcoll(9)
	character(len=50) :: filename
	logical, external :: readflags
	
	call chdir(rundir)
	!VN: For input readin
	open(file_input,file = "input.txt",status="old",access="sequential",action="read")

	if (readflags(file_input,'Peptide sequence 1')) then
		read(file_input,*) pep1
		write(fileout,*) 'Peptide sequence 1: ', pep1
	else
		write(fileout,*) 'Read input: FAIL TO READ PEPTIDE SEQUENCE 1'
		Stop
	endif
	if (readflags(file_input,'Number of residues in peptide 1')) then
		read(file_input,*) chnln1
		write(fileout,'(A,i7)') 'Number of residues in peptide 1: ', chnln1
	else
		write(fileout,*) 'Read input: FAIL TO READ NUMBER OF RESIDUES IN PEPTIDE 1'
		Stop
	endif
	if (readflags(file_input,'Number of Glycines in peptide 1')) then
		read(file_input,*) numgly1
		write(fileout,'(A,i2)') 'Number of Glycines in peptide 1: ', numgly1
	else
		write(fileout,*) 'Read input: FAIL TO READ NUMBER OF GLYCINES IN PEPTIDE 1'
		Stop
	endif
	if (readflags(file_input,'Number of peptide 1 chains in the system')) then
		read(file_input,*) nc
		write(fileout,'(A,i3)') 'Number of peptide 1 chains in the system: ', nc
	else
		write(fileout,*) 'Read input: FAIL TO READ NUMBER OF PEPTIDE 1 CHAINS IN THE SYSTEM'
		Stop
	endif
!!!!!!!!!!!!!!!!!
	if (readflags(file_input,'Peptide sequence 2')) then
		read(file_input,*) pep2
		write(fileout,*) 'Peptide sequence 2: ', pep2
	else
		write(fileout,*) 'Read input: FAIL TO READ PEPTIDE SEQUENCE 2'
		Stop
	endif
	if (readflags(file_input,'Number of residues in peptide 2')) then
		read(file_input,*) chnln2
		write(fileout,'(A,i2)') 'Number of residues in peptide 2: ', chnln2
	else
		write(fileout,*) 'Read input: FAIL TO READ NUMBER OF RESIDUES IN PEPTIDE 2'
		Stop
	endif
	if (readflags(file_input,'Number of Glycines in peptide 2')) then
		read(file_input,*) numgly2
		write(fileout,'(A,i2)') 'Number of Glycines in peptide 2: ', numgly2
	else
		write(fileout,*) 'Read input: FAIL TO READ NUMBER OF GLYCINES IN PEPTIDE 2'
		Stop
	endif
	if (readflags(file_input,'Number of peptide 2 chains in the system')) then
		read(file_input,*) nc2
		write(fileout,'(A,i3)') 'Number of peptide 2 chains in the system: ', nc2
	else
		write(fileout,*) 'Read input: FAIL TO READ NUMBER OF PEPTIDE 2 CHAINS IN THE SYSTEM'
		Stop
	endif
!!!!!!!!!!!!!!!!
	if (readflags(file_input,'Box length in Angstrom')) then
		read(file_input,*) boxlength
		write(fileout,'(A,f7.1)') 'Box length in Angstrom: ', boxlength
	else
		write(fileout,*) 'Read input: FAIL TO READ BOX LENGTH'
		Stop
	endif
	if (readflags(file_input,'Simulation temperature in Kelvin')) then
		read(file_input,*) simtemp
		write(fileout,'(A,f7.1)') 'Simulation temperature in Kelvin: ', simtemp
	else
		write(fileout,*) 'Read input: FAIL TO READ SIMULATION TEMPERATURE'
		Stop
	endif
	if (readflags(file_input,'Result recording frequency in collisions')) then
		read(file_input,*) simcoll
		write(fileout,'(A,i17)') 'Result recording frequency in collisions: ', simcoll
	else
		write(fileout,*) 'Read input: FAIL TO READ SIMULATION COLLISIONS'
		Stop
	endif
	!read(file_input,*) maxtemp
	!read(file_input,*) mintemp
	!read(file_input,*) annealcoll
	close(file_input)

	nb1=4*chnln1          		! number of beads
    	numbeads1=nb1-numgly1         ! number of beads without glycines
	nopwg1=nb1*nc      	! total number of beads
    	nop1=numbeads1*nc    ! total number of beads minus glycines
    	nb2=4*chnln2          		! number of beads
    	numbeads2=nb2-numgly2         ! number of beads without glycines
	nopwg2=nb2*nc2      	! total number of beads
    	nop2=numbeads2*nc2   ! total number of beads minus glycines
	
	redsimT = (simtemp+115.79)/2288.46 
	open(175,file='temperatures/simtemp',status='unknown')
	write(175,*) redsimT
	write(175,*) simcoll
	close(175)

	
	end subroutine readinputs
