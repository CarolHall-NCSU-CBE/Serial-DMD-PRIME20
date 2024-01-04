!	subroutine events finds type of event for i with each neighbor
!	j, sends pair to core, bond, or sqwell subroutine for 
!	time calculation, stores event time and partner for the
!	first event that will occur with each 

!	the point of this subr is to find tim, nptnr, and coltype for 
!	each nop
			   
!	j is always > i, but aa(j) not necessarily > aa(j)

	subroutine events()
 
#include "def.h"

	use global
	use inputreadin
	
	implicit none

	integer evcode,type,i,j,k,kstart,kend
	real*8 tij
	
	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do i=1,(noptotal)
		kstart=(i-1)*maxnbs+1             
		kend=kstart+na_npt(i)-1       
		do k=kstart,kend              
			j=nb(k)           
	      		tij=1000000000.d0
	      		evcode=ev_code(j,i)
	      		if (evcode.le.3) then
		 		call core(i,j,evcode,tij,type)
              	elseif ((evcode.ge.4).and.(evcode.le.12)) then
		 		call bond(i,j,evcode,tij,type)
	      		elseif (evcode.eq.15) then
		 		call nc_sqwel(i,j,evcode,tij,type)
	      		elseif (evcode.eq.16) then
		 		call sqwel(i,j,evcode,tij,type)
	      		elseif ((evcode.ge.17).and.(evcode.le.26))then
		 		call core(i,j,evcode,tij,type)
	      		elseif (evcode.ge.40) then
                 		call sqshlder(i,j,evcode,tij,type)
	      		else 
		 		write(6,*)'error in ev_code matrix (events)'
		 		write(6,*)'ev_code(i,j)=',evcode
		 		write(6,*)'i=',i
		 		write(6,*)'j=',j
		 		call exit(-1)
	      		endif

!             	at this point, have found event time for a given i,j pair.  
!             	test to see if it's the soonest event i will face; 
!             	if so, save that time, partner, and event type.
	      		if (tij.lt.tim(i)) then
		 		tim(i)=tij
		 		nptnr(i)=j
		 		coltype(i)=type 
	      		endif

#ifdef debugging
	      		if ((tij.lt.0.0) .and. (abs(tij) .gt. 1.0e-10)) then
		 		write(6,*)'in events'
		 		write(6,*)'coll=',coll
		 		write(6,*)'boxl=',boxl
		 		write(6,*)'factor=',ev_param(1,evcode)
		 		write(6,*)'i =',i
		 		write(6,*)'j =',j
		 		write(6,*)'tij=',tij
		 		write(6,*)'coltype=',type
		 		write(6,*)'evcode=',evcode
                 		write(6,*)'stop, tij is less than 0'
		 		call exit(-1)
	      		endif
#endif
		enddo
		do k = 1, 3
	      		j = extra_repuls(i,k)
	      		if (j .gt. i) then
		 		tij=1000000000.0
	         		evcode=ev_code(j,i)
	         		if (evcode == 1) then
	            			call core(i,j,evcode,tij,type)
	         		else
	            			call sqshlder(i,j,evcode,tij,type)
	         		end if
#ifdef debugging
				if ((tij.lt.0.0) .and. (abs(tij) .gt. 1.0e-10)) then
                    			write(6,*)'in events'
                    			write(6,*)'coll=',coll
                    			write(6,*)'boxl=',boxl
                    			write(6,*)'factor=',ev_param(1,evcode)
                    			write(6,*)'i =',i
                    			write(6,*)'j =',j
                    			write(6,*)'tij=',tij
                    			write(6,*)'coltype=',type   
                    			write(6,*)'stop, tij is less than 0 in events, extra_repulse', k
                    			call exit(-1)
                 		endif
#endif

                 		if (tij .lt. tim(i)) then
                    			nptnr(i)=j
                    			coltype(i)=type
                    			tim(i)=tij
                 		end if
              	endif
	   	enddo
	enddo
!write(6,*)i,j,chnnum(i),chnnum(j),identity(i),identity(j),tij
	do i = 1, numbin+1
		bin(i)=0
	end do

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do i = 1, (noptotal) + 3
		tlinks(i) = 0
		tlinks2(i) = 0
if (coll .eq. 2131749) then
write(6,*)i,nptnr(i)
write(6,*)identity(i),identity(nptnr(i))
write(6,*)tim(i)
endif
		if (tim(i) .lt. interval_max) call add_tbin(i)
	end do

	return

	end
	
