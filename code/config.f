	subroutine config()

#include "def.h"

	use global
	use inputreadin

	implicit none

    !LR: Changed a hardcoded 2-species variable reference to a noptotal variable	
	real*8 xtmp(noptotal),ytmp(noptotal),ztmp(noptotal)
	integer k
	character*64 filename

	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do k = 1,(noptotal)
		xtmp(k)=sv(1,k)+sv(4,k)*tfalse
		ytmp(k)=sv(2,k)+sv(5,k)*tfalse
		ztmp(k)=sv(3,k)+sv(6,k)*tfalse
	enddo

! VN: Write configuration file:
	write(runcf) coll, t+tfalse,xtmp, ytmp, ztmp
! VN: Write bond partner file: 
	write(runpartner) coll, bptnr

	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do k = 1,(noptotal)
		xtmp(k)=sv(4,k)
		ytmp(k)=sv(5,k)
		ztmp(k)=sv(6,k)
	enddo

! VN: Write velocity file:
	rewind(runlasvel)
	write(runlasvel) coll, xtmp, ytmp, ztmp

	
        return

        end
