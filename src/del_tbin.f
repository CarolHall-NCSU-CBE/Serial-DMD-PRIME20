	subroutine del_tbin(i)

#include "def.h"

	use global
	use inputreadin

	implicit none

	integer i

	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	if (tlinks2(i) .gt. (noptotal)+3) then
		bin(tlinks2(i)-(noptotal)-3) = tlinks(i)
	else
		tlinks(tlinks2(i))=tlinks(i)
	end if

	if (tlinks(i) .ne. 0) tlinks2(tlinks(i))=tlinks2(i)
	
	tlinks2(i) = 0

	return

	end
