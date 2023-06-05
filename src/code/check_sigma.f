	subroutine check_sigma(i,j,rating)

#include "def.h"

	use global
	use inputreadin

	implicit none
	
	integer i,j
	real*8 vxij,vyij,vzij,rxij,ryij,rzij,rijsq,rij,diff,rating

	vxij=sv(4,i)-sv(4,j)
      	vyij=sv(5,i)-sv(5,j)
      	vzij=sv(6,i)-sv(6,j)
      	rxij=sv(1,i)-sv(1,j) + vxij*tfalse
      	ryij=sv(2,i)-sv(2,j) + vyij*tfalse
      	rzij=sv(3,i)-sv(3,j) + vzij*tfalse
      	rxij=rxij-dnint(rxij)
      	ryij=ryij-dnint(ryij)
      	rzij=rzij-dnint(rzij)
      	rijsq=rxij*rxij+ryij*ryij+rzij*rzij
      	rij=dsqrt(rijsq)
      	diff=rijsq-sigma_sq(identity(i),identity(j))
      	if (diff.lt.0.d0) then
		rating = 10.d0
      	else
		rating = 11.d0
      	end if

      	return
      
      	end
