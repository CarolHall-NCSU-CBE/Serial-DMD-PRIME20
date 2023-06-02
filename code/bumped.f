      	subroutine bumpoff(i,j,evcode)
	
#include "def.h"

      	use global
	use inputreadin

      	implicit none
	
      	integer i,j,evcode
      	real vxij,vyij,vzij,rxij,ryij,rzij,bij,bumpdist

      	vxij=sv(4,i)-sv(4,j)
      	vyij=sv(5,i)-sv(5,j)
      	vzij=sv(6,i)-sv(6,j)
      	rxij=sv(1,i)-sv(1,j) + vxij*tfalse
      	ryij=sv(2,i)-sv(2,j) + vyij*tfalse
      	rzij=sv(3,i)-sv(3,j) + vzij*tfalse
      	rxij=rxij-dnint(rxij)
      	ryij=ryij-dnint(ryij)
      	rzij=rzij-dnint(rzij)
	bij=rxij*vxij+ryij*vyij+rzij*vzij
                     
	if (evcode .ge. 40) then
		bumpdist=smdist*dsqrt(shlddia_sq(identity(i),identity(j)))
	else
		bumpdist=smdist*dsqrt(welldia_sq(identity(i),identity(j)))
	end if

      	if (bij .lt. 0.d0) then
		sv(1,i) = sv(1,i) - bumpdist*rxij
           	sv(2,i) = sv(2,i) - bumpdist*ryij
           	sv(3,i) = sv(3,i) - bumpdist*rzij
           	sv(1,j) = sv(1,j) + bumpdist*rxij
           	sv(2,j) = sv(2,j) + bumpdist*ryij
           	sv(3,j) = sv(3,j) + bumpdist*rzij
	else
		sv(1,i) = sv(1,i) + bumpdist*rxij
           	sv(2,i) = sv(2,i) + bumpdist*ryij
           	sv(3,i) = sv(3,i) + bumpdist*rzij
		sv(1,j) = sv(1,j) - bumpdist*rxij
		sv(2,j) = sv(2,j) - bumpdist*ryij
		sv(3,j) = sv(3,j) - bumpdist*rzij
	endif

      	return
      
      	end
