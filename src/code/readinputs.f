#include "readflags.f90"
subroutine readinputs
	use inputreadin
	use global
	implicit none
	real annealtemp(9),  d1d2d3
	integer i, annealcoll(9),espos, annealcheck, error
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
	write(6,'(2A)') 'Peptide sequence 1: ', pep1
	
	chnln1 = len(trim(pep1))
	write(6,'(A,i3)') 'Number of residues in peptide 1: ', chnln1
	
	numgly1 = count(transfer(pep1,(/"x"/)).eq.'G')
	write(6,'(A,i3)') 'Number of Glycines in peptide 1: ', numgly1

	read(file_input,*) skipline	
	read(file_input,'(A)') readline
	espos = index(readline,'=')
	shortline = adjustl(readline((espos+1):len(readline)))
	read(shortline,*) nc
	write(6,'(A,i3)') 'Number of peptide 1 chains in the system: ', nc	

!!!!!!!!!!!!!!!!!
	read(file_input,*) skipline	
	read(file_input,'(A)') readline
	espos = index(readline,'=')
	pep2 = adjustl(readline((espos+1):len(readline)))
	write(6,'(2A)') 'Peptide sequence 2: ', pep2
	chnln2 = len(trim(pep2))
	write(6,'(A,i3)') 'Number of residues in peptide 2: ', chnln2
	numgly2 = count(transfer(pep2,(/"x"/)).eq.'G')
	write(6,'(A,i3)') 'Number of Glycines in peptide 2: ', numgly2

	read(file_input,*) skipline	
	read(file_input,'(A)') readline
	espos = index(readline,'=')
	shortline = adjustl(readline((espos+1):len(readline)))
	read(shortline,*) nc2
	write(6,'(A,i3)') 'Number of peptide 2 chains in the system: ', nc2	

!!!!!!!!!!!!!!!!
	read(file_input,*) skipline	
	read(file_input,'(A)') readline
	espos = index(readline,'=')
	shortline = adjustl(readline((espos+1):len(readline)))
	read(shortline,*) boxlength
	write(6,'(A,f7.1)') 'Box length in Angstrom: ', boxlength
	
	read(file_input,*) skipline	
	read(file_input,'(A)') readline
	espos = index(readline,'=')
	shortline = adjustl(readline((espos+1):len(readline)))
	read(shortline,*) simtemp
	write(6,'(A,f7.1)') 'Simulation temperature in Kelvin: ', simtemp
	simtemp = (simtemp+115.79)/2288.46 
	write(6,'(A,f7.3)') 'Simulation temperature in reduced unit: ', simtemp

	read(file_input,*) skipline	
	read(file_input,'(A)') readline
	espos = index(readline,'=')
	shortline = adjustl(readline((espos+1):len(readline)))
	read(shortline,*) simcoll
	write(6,'(A,i17)') 'Result recording frequency in collisions: ', simcoll
	
	read(file_input,*) skipline	
	read(file_input,'(A)') readline
	espos = index(readline,'=')
	shortline = adjustl(readline((espos+1):len(readline)))
	read(shortline,*) numsim

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
	read(shortline,*) newold

	if (newold == 0) then
		read(file_input,*) skipline	
		read(file_input,'(A)') readline
		espos = index(readline,'=')
		shortline = adjustl(readline((espos+1):len(readline)))
		read(shortline,*) annealcheck
		if (annealcheck == 0) then
			annealingsteps = 9
			allocate(temps(annealingsteps))
			temps= (/0.50,0.45,0.40,0.35,0.30,0.28,0.26,0.24,0.22/)	
			allocate(collset(annealingsteps))
			collset(1:5) = 100000000
			collset(6:7) = 150000000
			collset(8) = 200000000
			collset(9) = 250000000
		endif
		read(file_input,'(A)',iostat = error) readline
		if ((error == 0).and. (len(trim(readline)) .ne.0)) then
			espos = index(readline,'=')
			shortline = adjustl(readline((espos+1):len(readline)))
			read(shortline,*) d1d2d3
			dadjust1=d1d2d3
			dadjust2=d1d2d3
			write(6,'(A,f7.4)') 'Giving sidechain more space to find the correct location:', d1d2d3
			close(file_input)
		elseif ((error .ne. 0) .or.(len(trim(readline)) ==0) ) then  
			dadjust1=3.0
			dadjust2=3.0
			close(file_input)
		endif
	else
		close(file_input)
	endif
	end subroutine readinputs
