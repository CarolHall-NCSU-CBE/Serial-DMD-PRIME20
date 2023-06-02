	subroutine fileclose

#include "def.h"

	use global
	use inputreadin

	implicit none

!VN: Files that are stored in operating directory
	call chdir(rundir)
! Parameters:
	close(parahp1)
	close(parahp2)
	close(parafs1)
	close(parafs2)
	close(paraid) 
! Results:
	close(runcf)
	close(runpartner)
	close(runlasvel)
	close(rune)
	close(runphipsi)
	close(runpdb)	
	close(precf)
	close(runrca)
	close(prelasvel)
	close(fileout)

	
!VN: Files that are stored in source code package:
	call chdir(mydir)
	close(simwellha55a)
	close(simmass)
	close(simrtoall)
	close(simsqz)
	close(parabundles)
	close(pararn)
	close(pararc)
	close(parapro)
	close(parabdln)
	close(ep19ha55weak)
	
	return

	end
