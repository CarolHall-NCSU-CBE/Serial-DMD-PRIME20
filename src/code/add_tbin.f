
	subroutine add_tbin(i)

#include "def.h"

	use global
	use inputreadin

	implicit none

	integer i,j

	j = int((tim(i)+tbin_off)/sortsize)+1

#ifdef debugging

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	if (i .gt. (noptotal)+3) then
		print*, 'adding i>(nop1+nop2)+3', coll,i
	end if

	if (j .gt. numbin) then
		print*, j, ' is greater than numbin at', coll, i, coltype(i), tim(i)
		j = numbin
	end if

#endif

	tlinks(i)=bin(j)
	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	tlinks2(i)=j+(noptotal)+3
	if (bin(j) .ne. 0) tlinks2(bin(j)) = i
	bin(j)=i
	
	return

	end