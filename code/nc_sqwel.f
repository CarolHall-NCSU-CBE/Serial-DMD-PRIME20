!	subroutine nc_sqwel finds square well event time and type for 
!     	i,jpair; 4=sw capture, 5=sw dissociation, 6=sw bounce 7=break
!    	only 1,4 (or 4,1) pairs (unbound n and unbound c') 
!   	or 5,8 (or 8,5) pairs **in which they are bound to each other**
!   	may have square-well interaction, otherwise core collision

	subroutine nc_sqwel(i,j,evcode,tij,type)

#include "def.h"

     	use global
	use inputreadin

	implicit none

	integer evcode,type,i,j
	real*8 vxij,vyij,vzij,rxij,ryij,rzij,bij,rijsq,vijsq
     	real*8 diff,corediscr,fac_sigsq,fac_cored,welldiscr,tij

	vxij=sv(4,i)-sv(4,j)
	vyij=sv(5,i)-sv(5,j)	
	vzij=sv(6,i)-sv(6,j)
	rxij=sv(1,i)-sv(1,j)+vxij*tfalse
	ryij=sv(2,i)-sv(2,j)+vyij*tfalse
	rzij=sv(3,i)-sv(3,j)+vzij*tfalse
	rxij=rxij-dnint(rxij)
	ryij=ryij-dnint(ryij)
	rzij=rzij-dnint(rzij)
	bij=rxij*vxij+ryij*vyij+rzij*vzij
	rijsq=rxij*rxij+ryij*ryij+rzij*rzij
	vijsq=vxij*vxij+vyij*vyij+vzij*vzij
	diff=rijsq-welldia_sq(identity(i),identity(j))	
!   	first, consider unbound n and unbound c'

	if ((identity(i)+identity(j)).eq.5) then
		if (bij.lt.0.d0) then
!          		then i and j are moving together
	      		if (diff.lt.0.d0) then
!             		then wells are already overlapping
	         		corediscr=bij*bij-vijsq*(rijsq-sigma_sq(identity(i),identity(j)))
		 		if (corediscr.gt.0.d0) then
!                			then a core collision will occur
		    			tij=(-bij-dsqrt(corediscr))/vijsq
		    			type=1
		 		else
!                			then leave the well - remember to turn off the square shoulder potential
	            			welldiscr=bij*bij-vijsq*diff
                    			tij=(-bij+dsqrt(welldiscr))/vijsq
                    			type=16
		 		endif
	      		else
!             		the wells are not yet overlapping
	         		welldiscr=bij*bij-vijsq*diff
		 		if (welldiscr.gt.0.d0) then
!                			then a capture may occur (if the square shoulder particles are far apart)
!                			assume capture
		    			tij=(-bij-dsqrt(welldiscr))/vijsq
		    			type=7
!                			in main, right before this event, determine if it's a capture 
		 		endif
	      		endif
	   	else
!          		then leave the well - remember to turn off the square shoulder potential
              	if (diff.lt.0.d0) then
	         		welldiscr=bij*bij-vijsq*diff
                 		tij=(-bij+dsqrt(welldiscr))/vijsq
                 		type=16	 		
	      		endif
	   	endif
	elseif (bptnr(i).eq.j) then
!      	i and j are already bound to each other
	   	if (bij.lt.0.d0) then
!          		then i and j are moving together
	      		fac_sigsq=sigma_sq(identity(i),identity(j))*ev_param(1,15)*ev_param(1,15)
	      		fac_cored=bij*bij-vijsq*(rijsq-fac_sigsq)
	      		if (fac_cored.gt.0.d0) then
!             		then a core collision will occur
		 		tij=(-bij-dsqrt(fac_cored))/vijsq
		 		type=1
	      		else 
!                		cores are moving towards each other but will not touch
!                		either dissociation or bounce, times are the same
!                		in eventdyn, decide which event occurs
	         		welldiscr=bij*bij-vijsq*diff
		 		tij=(-bij+dsqrt(welldiscr))/vijsq
		 		type=8
	      		endif
	   	else
!          		i and j are moving away from each other
!          		either dissociation or bounce, times are the same
!          		in main, right before event, check angle and ho distance to
!          		decide which event occurs
	      		welldiscr=bij*bij-vijsq*diff
	      		tij=(-bij+dsqrt(welldiscr))/vijsq
	      		type=8

!if (coll .ge. 1500183) then
!write(fileout,*)
!write(fileout,*)'This is where my Negative collision occurs'
!write(fileout,*)'Collision:',coll
!write(fileout,*)'Time',tij
!write(fileout,*)'HB Partners:',i,bptnr(i)
!write(fileout,*)'I and J:',i,j
!write(fileout,*)'IDs:',identity(i),identity(j)
!write(fileout,*)'Chains:',chnnum(i),chnnum(j)
!write(fileout,*)'EV_CODE',ev_code(i,j)
!write(fileout,*)'Coltype',type
!endif
	   	endif
	else
!       	beads are not eligible for square-well interaction
!       	core collision instead
	   	if (bij.lt.0.d0) then
!          		then i and j are moving together
              	welldiscr=bij*bij-vijsq*diff
              	if (welldiscr.gt.0.d0) then
                 		tij=(-bij-dsqrt(welldiscr))/vijsq
                 		type=9
!                		in main, right before this event, determine if it's a capture
!                		or no event
              	endif
	   	endif
	endif

	return

	end
