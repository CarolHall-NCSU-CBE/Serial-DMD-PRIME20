! =========================================================================
! This subroutine read inputs from text file
! This is a serial job
! Last modified on 10/12/2023 by Van Nguyen
! =========================================================================
#include "readflags.f90"
subroutine readinputs
	use inputreadin
	use global
	implicit none
	real d1d2d3
	integer i,espos, error, j
	character(len=50) :: filename,shortline
	character(len=175) :: skipline, readline
	logical, external :: readflags
	
	call chdir(rundir)
	!VN: For input readin
	open(file_input,file = "input.txt",status="old",access="sequential",action="read")

	read(file_input,*) skipline	
	read(file_input,'(A)') readline
	espos = index(readline,'=')
	pep1 = adjustl(readline((espos+1):len(readline)))
	write(6,'(2A)') ' Peptide sequence 1: ', pep1
	
	chnln1 = len(trim(pep1))
	write(6,'(A,i3)') ' Number of residues in peptide 1: ', chnln1
	
	numgly1 = count(transfer(pep1,(/"x"/)).eq.'G')
	write(6,'(A,i3)') ' Number of Glycines in peptide 1: ', numgly1

	read(file_input,*) skipline	
	read(file_input,'(A)') readline
	espos = index(readline,'=')
	shortline = adjustl(readline((espos+1):len(readline)))
	read(shortline,*) nc
	write(6,'(A,i3)') ' Number of peptide 1 chains in the system: ', nc	

!!!!!!!!!!!!!!!!!
	read(file_input,*) skipline	
	read(file_input,'(A)') readline
	espos = index(readline,'=')
	pep2 = adjustl(readline((espos+1):len(readline)))
	write(6,'(2A)') ' Peptide sequence 2: ', pep2
	chnln2 = len(trim(pep2))
	write(6,'(A,i3)') ' Number of residues in peptide 2: ', chnln2
	numgly2 = count(transfer(pep2,(/"x"/)).eq.'G')
	write(6,'(A,i3)') ' Number of Glycines in peptide 2: ', numgly2

	read(file_input,*) skipline	
	read(file_input,'(A)') readline
	espos = index(readline,'=')
	shortline = adjustl(readline((espos+1):len(readline)))
	read(shortline,*) nc2
	write(6,'(A,i3)') ' Number of peptide 2 chains in the system: ', nc2	

!!!!!!!!!!!!!!!!
	read(file_input,*) skipline	
	read(file_input,'(A)') readline
	espos = index(readline,'=')
	shortline = adjustl(readline((espos+1):len(readline)))
	read(shortline,*) boxlength
	write(6,'(A,f7.1)') ' Box length in Angstrom: ', boxlength
	
	read(file_input,*) skipline	
	read(file_input,'(A)') readline
	espos = index(readline,'=')
	shortline = adjustl(readline((espos+1):len(readline)))
	read(shortline,*) simtemp
	write(6,'(A,f7.1)') ' Simulation temperature in Kelvin: ', simtemp
	simtemp = (simtemp+115.79)/2288.467 
	write(6,'(A,f7.3)') ' Simulation temperature in reduced unit: ', simtemp

	read(file_input,*) skipline	
	read(file_input,'(A)') readline
	espos = index(readline,'=')
	shortline = adjustl(readline((espos+1):len(readline)))
	read(shortline,*) simcoll
	write(6,'(A,i17)') ' Result recording frequency in collisions: ', simcoll
	
	read(file_input,*) skipline	
	read(file_input,'(A)') readline
	espos = index(readline,'=')
	shortline = adjustl(readline((espos+1):len(readline)))
	read(shortline,*) trajcoll
	write(6,'(A,i17)') ' Result recording frequency in collisions: ', trajcoll
	
	nb1=4*chnln1          		! number of beads
    	numbeads1=nb1-numgly1         ! number of beads without glycines
	nopwg1=nb1*nc      	! total number of beads
    	nop1=numbeads1*nc    ! total number of beads minus glycines
    	nb2=4*chnln2          		! number of beads
    	numbeads2=nb2-numgly2         ! number of beads without glycines
	nopwg2=nb2*nc2      	! total number of beads
    	nop2=numbeads2*nc2   ! total number of beads minus glycines

		read(file_input,*) skipline	
		read(file_input,'(A)') readline
		espos = index(readline,'=')
		shortline = adjustl(readline((espos+1):len(readline)))
		read(shortline,*) annealcheck
		if (annealcheck == 1) then
			write(6,'(A)') ' Annealing process is run with temperature range defined by user.'
			read(file_input,'(A)',iostat = error) readline
			if ((error == 0) .and. (len(trim(readline)) .ne.0)) then
				espos = index(readline,'=')
				shortline = adjustl(readline((espos+1):len(readline)))
				read(shortline,*) maxannealing
				write(6,'(A,f7.1)') ' Starting of annealing temperature in Kelvin: ',  maxannealing
			endif
			read(file_input,'(A)',iostat = error) readline
			if ((error == 0) .and. (len(trim(readline)) .ne.0)) then
				espos = index(readline,'=')
				shortline = adjustl(readline((espos+1):len(readline)))
				read(shortline,*) minannealing
				write(6,'(A,f7.1)') ' Ending of annealing temperture in Kelvin: ',  minannealing
			endif
			read(file_input,'(A)',iostat = error) readline
			if ((error == 0) .and. (len(trim(readline)) .ne.0)) then
				espos = index(readline,'=')
				shortline = adjustl(readline((espos+1):len(readline)))
				read(shortline,*) incrannealing
				write(6,'(A,f7.1)') ' Temperature step for annealing in Kelvin: ',  incrannealing
			endif
			read(file_input,'(A)',iostat = error) readline
			if ((error == 0) .and. (len(trim(readline)) .ne.0)) then
				espos = index(readline,'=')
				shortline = adjustl(readline((espos+1):len(readline)))
				read(shortline,*) annealcoll
				write(6,'(A,i17)') ' Annealing collision: ',  annealcoll
			endif
		elseif (annealcheck == 0) then
			write(6,'(A)') ' Annealing process is run with default temperature from 1028K to 388K.'
		endif
		read(file_input,'(A)',iostat = error) readline
		if ((error == 0).and. (len(trim(readline)) .ne.0)) then
			espos = index(readline,'=')
			shortline = adjustl(readline((espos+1):len(readline)))
			read(shortline,*) d1d2d3
			dadjust1=d1d2d3
			dadjust2=d1d2d3
			write(6,'(A,f7.4)') ' Giving sidechain more space to find the correct location:', d1d2d3
			close(file_input)
		elseif ((error .ne. 0) .or.(len(trim(readline)) ==0) ) then  
			dadjust1=3.0
			dadjust2=3.0
			close(file_input)
		endif
		close(file_input)
	end subroutine readinputs
