!	subroutine eventredo takes (m,i) or (m,j) pair from main 
!	program and calculates the collision time for that pair. 
!	i and j were just involved in the collision. 
!	"i" appears on "m"'s neighbor list but "m" is not on "i"'s list.
!	"j" appears on "m"'s neighbor list but "m" is not on "j"'s list.
!	(due of the effort to avoid redundancy on the neighbor lists,
!	neighbor lists only include particles with higher numbers.
!	this subroutine is for neighbors with lower numbers.)
!
!	the point of this subr is to find tim, nptnr, and coltype for 
!	each nop
!	j is always > i, but aa(j) not necessarily > aa(j)

	subroutine eventredo_down(i,j) 

#include "def.h"

	use global
	use inputreadin

	implicit none

	integer evcode,type,i,j
	real*8 tij

	tij=1000000000.d0
	evcode=ev_code(i,j)

	if (evcode.le.3) then
		call core(i,j,evcode,tij,type)
	elseif ((evcode.ge.4).and.(evcode.le.12)) then
		call bond(i,j,evcode,tij,type)
!		note, bdln and cbl were scaled in scale_down.f
	elseif (evcode.eq.15) then
	   	call nc_sqwel(i,j,evcode,tij,type)
	elseif (evcode.eq.16) then
	   	call sqwel(i,j,evcode,tij,type)
	elseif ((evcode.ge.17).and.(evcode.le.26)) then
	   	call core(i,j,evcode,tij,type)
	elseif (evcode.ge.40) then
           	call sqshlder(i,j,evcode,tij,type)
	else 
	   	write(fileout,*)'error in ev_code matrix (events)'
	   	write(fileout,*)'ev_code(i,j)=',evcode
	   	write(fileout,*)'i=',i
	   	write(fileout,*)'j=',j
	endif

!       at this point, have found event time for a given i,j pair.  
!       test to see if it's the soonest event i will face; 
!       if so, save that time, partner, and event type.

#ifdef debugging

	if ((tij.lt.0.0) .and. (abs(tij) .gt. 1.0e-10)) then
	   	write(fileout,*)'in eventredo_down'
	   	write(fileout,*)'coll=',coll
	   	write(fileout,*)'boxl=',boxl
	   	write(fileout,*)'factor=',ev_param(1,evcode)
	   	write(fileout,*)'i =',i
	   	write(fileout,*)'j =',j
	   	write(fileout,*)'tij=',tij
	   	write(fileout,*)'coltype=',type
	  	write(fileout,*)'ev_code=',evcode
	   	write(fileout,*)'stop in eventredo_down, tij is less than 0'
	   	call exit(-1)
	endif

#endif

	tij = tij + tfalse

	if (tij.lt.tim(i)) then
		if (tlinks2(i) .ne. 0) call del_tbin(i)
	   	tim(i)=tij
	   	nptnr(i)=j
	   	coltype(i)=type
	   	if (tim(i) .lt. interval_max) call add_tbin(i)
	endif

	return

	end
	
