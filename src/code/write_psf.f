        subroutine write_psf

#include "def.h"

        USE GLOBAL
	use inputreadin
                      
	IMPLICIT NONE

        integer :: oxygen, i, j,izero, massid(23)
	character*1:: skipline,chainid
	character*4:: atom,segment(noptotal+chnln1*nc+chnln2*nc2), an(noptotal+chnln1*nc+chnln2*nc2)
	character*3:: resname(noptotal+chnln1*nc+chnln2*nc2)
	character*2 :: massresname(23)
	integer :: asn(noptotal+chnln1*nc+chnln2*nc2),resseqnum(noptotal+chnln1*nc+chnln2*nc2)
	real :: reducedmass(23), mass(23), masspsf(noptotal+chnln1*nc+chnln2*nc2)
	call chdir(mydir)
	open(unit=627, file = 'parameters/mass.data', status = 'unknown')
	do i = 1,23
		read(627,*) massresname(i), massid(i), reducedmass(i), mass(i)
	enddo
	call chdir(rundir)
	open(unit=773,file='prime20.psf',status='unknown')     
7       format(A4,3X,I4,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3)
17	format(A4,3X,I4,1X,A4,1X,A3,1X,A1,I4)
722	format(A6,3F9.3,3F7.2,X,A11,I4)
777	format(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,2G14.6)
	open(unit=732,file='results/run0001.pdb',status='old')
	read(732,*) skipline
	do i = 1, (noptotal+chnln1*nc+chnln2*nc2)
		read(732,17) atom,asn(i),an(i),resname(i),chainid,resseqnum(i)
	enddo
       	oxygen = chnln1*nc+chnln2*nc2
        write(773,'(A3)') 'PSF'
	write(773,*) ' '
	write(773, '(I8,1X,A6)') noptotal+oxygen,'!NATOM'
	izero = ichar('0')
	do j = 1,nc
		do i = 1,(chnln1*5-numgly1)
			segment((j-1)*(chnln1*5-numgly1)+i) = char(j+izero)//'P1'
		enddo
	enddo
	do j = 1,nc2
		do i = 1,(chnln2*5-numgly2)
			segment(nc*(chnln1*5-numgly1)+(j-1)*(chnln1*5-numgly1)+i) = char(j+izero)//'P2'
		enddo
	enddo
	do i = 1,(noptotal+oxygen)
		if (adjustl(an(i)) .eq. 'N') then
			masspsf(i) = mass(2)
		elseif (adjustl(an(i)) .eq. 'CA') then
			masspsf(i) = mass(1)
		elseif (adjustl(an(i)) .eq. 'C') then
			masspsf(i) = 12.011
		elseif (adjustl(an(i)) .eq. 'O') then
			masspsf(i) = mass(3)-12.011
		elseif (adjustl(an(i)) .eq. 'CB') then
			if (adjustl(resname(i)) .eq. 'GLY') then
				masspsf(i) = mass(4)
			elseif (adjustl(resname(i)) .eq. 'ARG') then
				masspsf(i) = mass(5)
			elseif (adjustl(resname(i)) .eq. 'ASN') then
				masspsf(i) = mass(6)
			elseif (adjustl(resname(i)) .eq. 'ASP') then
				masspsf(i) = mass(7)
			elseif (adjustl(resname(i)) .eq. 'GLN') then
				masspsf(i) = mass(8)
			elseif (adjustl(resname(i)) .eq. 'GLU') then
				masspsf(i) = mass(9)
			elseif (adjustl(resname(i)) .eq. 'HIS') then
				masspsf(i) = mass(10)
			elseif (adjustl(resname(i)) .eq. 'LYS') then
				masspsf(i) = mass(11)
			elseif (adjustl(resname(i)) .eq. 'PRO') then
				masspsf(i) = mass(12)
			elseif (adjustl(resname(i)) .eq. 'SER') then
				masspsf(i) = mass(13)
			elseif (adjustl(resname(i)) .eq. 'THR') then
				masspsf(i) = mass(14)
			elseif (adjustl(resname(i)) .eq. 'ALA') then
				masspsf(i) = mass(15)
			elseif (adjustl(resname(i)) .eq. 'CYS') then
				masspsf(i) = mass(16)
			elseif (adjustl(resname(i)) .eq. 'ILE') then
				masspsf(i) = mass(17)
			elseif (adjustl(resname(i)) .eq. 'LEU') then
				masspsf(i) = mass(18)
			elseif (adjustl(resname(i)) .eq. 'MET') then
				masspsf(i) = mass(19)
			elseif (adjustl(resname(i)) .eq. 'PHE') then
				masspsf(i) = mass(20)
			elseif (adjustl(resname(i)) .eq. 'TRP') then
				masspsf(i) = mass(21)
			elseif (adjustl(resname(i)) .eq. 'TYR') then
				masspsf(i) = mass(22)
			elseif (adjustl(resname(i)) .eq. 'VAL') then
				masspsf(i) = mass(23)
			endif
		endif
	enddo
	do i = 1,(noptotal+chnln1*nc+chnln2*nc2)
		write(773,'(I8,1X,A4,I4,1X,A4,1X,2A4,f14.6,f13.6,I7)') i, segment(i),resseqnum(i),resname(i),an(i),an(i),0.0,masspsf(i),0
	enddo
        end		