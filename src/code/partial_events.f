!	subroutine partial_events updates times for collisions btwn 
!    	pairs of  neighbors that include at least one of the two 
!    	particles that just collided (i and/or j) 

	subroutine partial_events(i,j,xpulse_del)

#include "def.h"

   	use global
	use inputreadin

	implicit none
	
	logical xpulse_del
	integer i,j,k,kk,l,ll,m,kstart,kend

	if (tlinks2(i) .ne. 0) call del_tbin(i)
    	tim(i)=interval_max+ltstep-tfalse
     	coltype(i)=-1         
	nptnr(i)=-1	
      	kstart=npt(i)
      	kend=kstart+na_npt(i)-1

      	do k=kstart,kend
     		l=nb(k)	
     		call eventredo_up(i,l)
	enddo

	do k = 1, 3
      		if (extra_repuls(i,k) .gt. i) then
         		call eventredo_up(i,extra_repuls(i,k))
      		endif
	end do

	tim(i)=tim(i)+tfalse
	if (tim(i) .lt. interval_max) call add_tbin(i)	

#ifdef canon
	if (j.ne.0) then
#endif	 
!   	recalculate for j's neighbors
	if (tlinks2(j) .ne. 0) call del_tbin(j)
    	tim(j)=interval_max+ltstep-tfalse
     	coltype(j)=-1
	nptnr(j)=-1
       kstart=(j-1)*maxnbs+1
  	kend=kstart+na_npt(j)-1

      	do k=kstart,kend
    		l=nb(k)
      		call eventredo_up(j,l)
    	enddo

	do k = 1, 3
     		if (extra_repuls(j,k) .gt. j) then
         		call eventredo_up(j,extra_repuls(j,k))
        	end if
	end do

	tim(j)=tim(j)+tfalse
	if (tim(j) .lt. interval_max) call add_tbin(j)	

#ifdef canon
	endif
#endif	

!   	recalculate for i's down neighbor list
!   	since i and j redirect after the event, these pairings are no longer valid
	kstart=(i-1)*maxnbs+1
	kend=kstart+nnabdn(i)-1
	
	do kk=kstart,kend 
		l=dnnab(kk)
!     		l's next event is *not* with i
		if (nptnr(l).ne.i) then
	   		call eventredo_down(l,i)
!     			l's next event is with i
		else
	      		if (tlinks2(l) .ne. 0) call del_tbin(l)
	      		tim(l)=interval_max+ltstep-tfalse
	      		coltype(l)=-1
	      		nptnr(l)=-1
              	kstart=(l-1)*maxnbs+1
              	kend=kstart+na_npt(l)-1
              	do ll=kstart,kend
                 		m=nb(ll)
                 		call eventredo_up(l,m)
              	enddo
	      		do ll=1,3
                 		if (extra_repuls(l,ll) .gt. l) then
	            			call eventredo_up(l,extra_repuls(l,ll))
                 		end if
	      		end do
	      		tim(l)=tim(l)+tfalse
	      		if (tim(l) .lt. interval_max) call add_tbin(l)
	   	endif
	enddo

	do kk = 1, 3
		l = extra_repuls(i,kk)
		if ((l .lt. i) .and. (l .ne. 0)) then
	      		if (nptnr(l).ne.i)  then
	         		call eventredo_down(l,i)
	      		else
	        		if (tlinks2(l) .ne. 0) call del_tbin(l)
                 		tim(l)=interval_max+ltstep-tfalse
	         		coltype(l)=-1
                 		nptnr(l)=-1
                 		kstart=(l-1)*maxnbs+1             
                 		kend=kstart+na_npt(l)-1       
                 		do ll=kstart,kend              
                    			m=nb(ll)           
                    			call eventredo_up(l,m)             
                 		enddo         
                 		do ll = 1, 3
                    			if (extra_repuls(l,ll) .gt. l) then
                       			call eventredo_up(l,extra_repuls(l,ll))
                    			end if
                 		end do
                 		tim(l)=tim(l)+tfalse
		 		if (tim(l) .lt. interval_max) call add_tbin(l)
	      		endif	 
           	end if
	end do

#ifdef canon
	if (j.ne.0) then
#endif

!    	recalculate for j's down neighbor list
	kstart=(j-1)*maxnbs+1
	kend=kstart+nnabdn(j)-1

	do kk=kstart,kend
		l=dnnab(kk)
!       	l is not i 
	   	if (l.ne.i) then
!             	and l's next event is *not* with j
	      		if (nptnr(l).ne.j)  then
		 		call eventredo_down(l,j)
!             		(l is still not i) but  l's next event is with j
	      		else
	         		if (tlinks2(l) .ne. 0) call del_tbin(l)
		 		tim(l)=interval_max+ltstep-tfalse
		 		coltype(l)=-1
                 		nptnr(l)=-1
                 		kstart=(l-1)*maxnbs+1             
                 		kend=kstart+na_npt(l)-1       
                 		do ll=kstart,kend              
                    			m=nb(ll)           
                    			call eventredo_up(l,m)          
                 		enddo         
	         		do ll = 1, 3
                    			if (extra_repuls(l,ll) .gt. l) then
	               			call eventredo_up(l,extra_repuls(l,ll)) 
                    			end if
	         		end do

                 		tim(l)=tim(l)+tfalse
	         		if (tim(l) .lt. interval_max) call add_tbin(l)
	      		endif
	   	endif
	enddo

    	do kk = 1, 3  
           	l = extra_repuls(j,kk)
           	if ((l .lt. j) .and. (l .ne. 0)) then
              	if (nptnr(l).ne.j)  then
                 		call eventredo_down(l,j)
              	else
	         		if (tlinks2(l) .ne. 0) call del_tbin(l)
                 		tim(l)=interval_max+ltstep-tfalse
                 		coltype(l)=-1
                 		nptnr(l)=-1
	         		kstart=(l-1)*maxnbs+1             
                 		kend=kstart+na_npt(l)-1       
                 		do ll=kstart,kend              
                    			m=nb(ll)           
                    			call eventredo_up(l,m)             
                 		enddo         
                 		do ll = 1, 3
                    			if (extra_repuls(l,ll) .gt. l) then
                       			call eventredo_up(l,extra_repuls(l,ll))
                    			end if
                 		end do
                 		tim(l)=tim(l)+tfalse
                 		if (tim(l) .lt. interval_max) call add_tbin(l)
              	endif 
           	end if
        end do

#ifdef canon
	endif
#endif 

	if (xpulse_del) then
	   	if (identity(i) .lt. identity(j)) then
	      		call repuls_del_b(i,j)
	   	else
	      		call repuls_del_b(j,i)
  	   	endif
	endif

	return

	end
