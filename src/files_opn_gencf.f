! VN: This subroutine is to open existing files or creating new files needed for initial configuration generating step
	subroutine gencf_file_opener

#include "def.h"

	use global
	use inputreadin

	implicit none

	integer i, izero
	character*64 filename, zfilename
	logical exist_flag, zexist_flag
	character*1 zero
	character*1 char
	integer ichar, len
	data zero/'0'/

!VN: Files that are stored in operating directory
	call chdir(rundir)
!VN: For input readin
	open(file_input,file = "input.txt",status="old",access="sequential",action="read")
! First peptide sequence ID:
	open(ideninp,file='inputs/identity.inp',status='unknown')
! Second peptide sequence ID:
	open(iden2inp,file='inputs/identity2.inp',status='unknown')
! Check:
	open(checkid,file='checks/identity.out')
	open(checkcf1,file='checks/configone1.pdb')
	open(checkcf2,file='checks/configone2.pdb')
	open(checkcf,file='checks/config.pdb')
	open(checkcfout,file='checks/configcheck.out')
	open(checkmass,file='checks/masses.out')
	open(checksumvel,file='checks/sumvelcheck.out')
	open(checkvel,file='checks/velocitycheck.out')
	open(checkhb,file='checks/hbcheck.out')



! Parameters:
	open(parahp1,file='parameters/hp1.inp')
	open(parahp2,file='parameters/hp2.inp')
	open(parafs1,file='parameters/firstside1.data')
	open(parafs2,file='parameters/firstside2.data')
	open(paraid,file='parameters/identity.inp') 
! Initial configurations:
	open(inpinfon1,file='inputs/chninfo-n1.data',status='unknown',form='formatted')
	open(inpinfon2,file='inputs/chninfo-n2.data',status='unknown',form='formatted')

! Results:
	izero = ichar('0')
	i = 0
	
	!filename = 'results/run'//'0000.energy'
	!zfilename = 'results/run'//'0000.energy'//'.gz'
	fname_digits = '0000'
	open(runcf,file='results/run'//fname_digits//'.config',status='unknown',form='unformatted')
	open(runlasvel,file='results/run'//fname_digits//'.lastvel',status='unknown',form='unformatted')
	open(runpartner,file='results/run'//fname_digits//'.bptnr',status='unknown',form='unformatted')
	open(rune,file='results/run'//fname_digits//'.energy',status='unknown',form='unformatted')

!VN: Files that are stored in source code package:
	call chdir(mydir)
	open(genericpepx,file='genconfig/inputs/peptidex.inp',status='old')  
	open(genericpepy,file='genconfig/inputs/peptidey.inp',status='old')  
	open(genericpepz,file='genconfig/inputs/peptidez.inp',status='old')
	open(ffha55a,file='genconfig/inputs/beadwell_ha55a.data',status='unknown')
	open(massfile,file='genconfig/parameters/mass.data',status = 'old')
	open(parartoall,file='genconfig/parameters/rcarnrco.data',status='old')
	OPEN(parasqz,file='genconfig/parameters/sqz6to10.data',STATUS='OLD')
	
	return

	end

