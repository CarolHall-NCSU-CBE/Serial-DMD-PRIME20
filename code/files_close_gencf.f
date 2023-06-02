	subroutine gencf_filecloser

#include "def.h"

	use global
	use inputreadin

	implicit none

!VN: Files that are stored in operating directory
	call chdir(rundir)
! First peptide sequence ID:
	close(ideninp)
! Second peptide sequence ID:
	close(iden2inp)
! Check:
	close(checkid)
	close(checkcf1)
	close(checkcf2)
	close(checkcf)
	close(checkcfout)
	close(checkmass)
	close(checksumvel)
	close(checkvel)
	close(checkhb)

! Parameters:
	close(parahp1)
	close(parahp2)
	close(parafs1)
	close(parafs2)
	close(paraid) 
! Initial configurations:
	close(inpinfon1)
	close(inpinfon2)

! Results:
	close(runcf)
	close(runlasvel)
	close(runpartner)
	close(rune)

!VN: Files that are stored in source code package:
	call chdir(mydir)
	close(genericpepx)  
	close(genericpepy) 
	close(genericpepz)
	close(ffha55a)
	close(massfile)
	close(parartoall)
	close(parasqz)
	
	return

	end
