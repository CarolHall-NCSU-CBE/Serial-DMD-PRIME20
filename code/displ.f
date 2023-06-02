!	subroutine displ tracks mc displacement between neighbor list 
!	updates. also, after each collision it determines whether the 
!	list should be updated.
 
	subroutine displ(update)

#include "def.h"

	use global
	use inputreadin

	implicit none

	logical update
	real*8 dis, moved_far,a,b,c,moved
	integer i,j,k
		
!	determines the vector from the last position to the current 
!	position for each molecule 

	moved_far = 0.d0

	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do i=1,(noptotal)
		a=old_rx(i)-sv(1,i)
		b=old_ry(i)-sv(2,i)
		c=old_rz(i)-sv(3,i)
!	   	find total displacement since last neighbor list update
	   	dis = a*a+b*b+c*c
	   	moved=dis/hdelr
	   	if (moved .gt. moved_far) then 
	      		moved_far = moved
	   	endif
	end do

!	detm if neighbor list should be updated after current collision

	if (moved_far .ge. 0.1d0) then
		update = .true.
		if (moved_far.ge.1.25d0**2) then
			t_fact = t_fact / 1.01d0
			interval = t_fact/dsqrt(setemp)
			print*, 'moved too far by', moved_far, 't_fact', t_fact, ' at t', t+tfalse
		endif
	else
		update=.false.
	endif

	return

	end
	
	
	
