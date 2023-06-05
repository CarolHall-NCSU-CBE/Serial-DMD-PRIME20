!	subroutine eventdyn calculates new velocities for particles 
!	involved in a collision and calculates collision virial
!	the notations in PRIME20 and the assigment are different

	subroutine eventdyn (i,j,evcode) 

#include "def.h"

	use global
	use inputreadin

	implicit none

	integer evcode,i,j,k,ii,jj,ie
	real*8 vxij,vyij,vzij,rxij,ryij,rzij,bij,sigsq,ratio,blmin,blmax,del_pe
	real*8 wellsq,epsave,delvx,delvy,delvz,bumpdist
	real*8 rmass

!	apply periodic boundary conditions to i,j pair
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
	rmass=2*bm(i)*bm(j)/(bm(i)+bm(j))
	if(rmass.le.0) write(fileout,*) "rmass=",rmass
	if(rmass.le.0) write(fileout,*)  i,j,identity(i),identity(j)
	if(bm(i).le.0) write(fileout,*) "bm(i)=",bm(i)
	if(rmass.le.0) write(fileout,*)  i,identity(i)
	if(bm(j).le.0) write(fileout,*) "bm(j)=",bm(j)
	if(rmass.le.0) write(fileout,*)  j,identity(j)

       if (i .le. nop1) then
		if (coltype(i).eq.2) then
!       		bond collision
	   		blmin=ev_param(2,evcode)
			ii=i-((chnnum(i)-1)*numbeads1)
	   		if (evcode.eq.10) then
              		blmin=(1.d0-del_bdln(ii))*bdln(ii)
	   		elseif (evcode.eq.11) then
              		blmin=(1.d0-del_blrn(ii-chnln1))*bl_rn(ii-chnln1)
	   		elseif (evcode.eq.12) then		
              		blmin=(1.d0-del_blrc(ii-2*chnln1))*bl_rc(ii-2*chnln1)
	   		endif
			ratio=rmass*bij/(blmin*blmin)	   
		elseif (coltype(i).eq.3) then
!       		bond extension
	   		blmax=ev_param(3,evcode)
	   		ii=i-((chnnum(i)-1)*numbeads1)
	   		if (evcode.eq.10) then
              		blmax=(1.d0+del_bdln(ii))*bdln(ii)
	   		elseif (evcode.eq.11) then
              		blmax=(1.d0+del_blrn(ii-chnln1))*bl_rn(ii-chnln1)
	   		elseif (evcode.eq.12) then		
              		blmax=(1.d0+del_blrc(ii-2*chnln1))*bl_rc(ii-2*chnln1)
	   		endif
	   		ratio=rmass*bij/(blmax*blmax)
		elseif (coltype(i).eq.1) then
!       		core collision	
           		if (evcode .eq. 15) then
              		if (bptnr(i) .eq. j) then
                 			sigsq = sigma_sq(identity(i),identity(j))*ev_param(1,evcode)**2
              		else
                 			sigsq = sigma_sq(identity(i),identity(j))
              		end if
           		else
              		sigsq = sigma_sq(identity(i),identity(j))*ev_param(1,evcode)**2
              		if(evcode.ge.22.and.evcode.le.26) then
              			k=max0(identity(i),identity(j))
              			sigsq=sigsq*sqz610(evcode-21,k)**2
              		endif
           		end if
	   		ratio=rmass*bij/sigsq
        	else if (coltype(i).eq.4) then
	   		wellsq=welldia_sq(identity(i),identity(j))
           		epsave=ep_sqrt(identity(i),identity(j))
           		del_pe = 4.0*wellsq*epsave/rmass
		       if (bij*bij+del_pe.gt.0.d0) then
!          			capture
	   			ratio=rmass*(dsqrt((4.d0*wellsq*epsave/rmass)+bij*bij)+bij)/(2.d0*wellsq)
	   			coltype(i) = 20
	   			bumpdist=smdist*dsqrt(wellsq)
           			sv(1,i) = sv(1,i) - bumpdist*rxij
           			sv(2,i) = sv(2,i) - bumpdist*ryij
           			sv(3,i) = sv(3,i) - bumpdist*rzij
           			sv(1,j) = sv(1,j) + bumpdist*rxij
           			sv(2,j) = sv(2,j) + bumpdist*ryij
           			sv(3,j) = sv(3,j) + bumpdist*rzij
           		else    
! 				bounce
	      			ratio=rmass*bij/wellsq
              		coltype(i)=22
	   			bumpdist=smdist*dsqrt(wellsq)
           			sv(1,i) = sv(1,i) + bumpdist*rxij
           			sv(2,i) = sv(2,i) + bumpdist*ryij
           			sv(3,i) = sv(3,i) + bumpdist*rzij
           			sv(1,j) = sv(1,j) - bumpdist*rxij
           			sv(2,j) = sv(2,j) - bumpdist*ryij
           			sv(3,j) = sv(3,j) - bumpdist*rzij
           		endif
		else if (coltype(i)==8) then
                  	wellsq=welldia_sq(identity(i),identity(j))
           		epsave=ep_sqrt(identity(i),identity(j))
           		del_pe = 4.0*wellsq*epsave/rmass
	   		bumpdist=smdist*dsqrt(wellsq)
           		if (bij*bij.gt.del_pe) then
!          			k.e. high enough for dissociation event
	      			ratio=rmass*(-dsqrt(-del_pe+bij*bij)+bij)/(2.d0*wellsq)
              		coltype(i)=21		       
              		sv(1,i) = sv(1,i) + bumpdist*rxij
              		sv(2,i) = sv(2,i) + bumpdist*ryij
              		sv(3,i) = sv(3,i) + bumpdist*rzij
              		sv(1,j) = sv(1,j) - bumpdist*rxij
              		sv(2,j) = sv(2,j) - bumpdist*ryij  
              		sv(3,j) = sv(3,j) - bumpdist*rzij
!          			k.e. not high enough, therefore they bounce
           		else
	      			ratio=rmass*bij/wellsq
              		coltype(i)=22
              		sv(1,i) = sv(1,i) - bumpdist*rxij
              		sv(2,i) = sv(2,i) - bumpdist*ryij
              		sv(3,i) = sv(3,i) - bumpdist*rzij
              		sv(1,j) = sv(1,j) + bumpdist*rxij
              		sv(2,j) = sv(2,j) + bumpdist*ryij
              		sv(3,j) = sv(3,j) + bumpdist*rzij
           		endif
!       		core collision at square-well diameter due to overcrowding 
		else if (coltype(i)==9) then 
	   		wellsq=welldia_sq(identity(i),identity(j))
	   		ratio=rmass*bij/wellsq
           		coltype(i)=23
	   		bumpdist=smdist*dsqrt(wellsq)
           		sv(1,i) = sv(1,i) + bumpdist*rxij
           		sv(2,i) = sv(2,i) + bumpdist*ryij
           		sv(3,i) = sv(3,i) + bumpdist*rzij
           		sv(1,j) = sv(1,j) - bumpdist*rxij
           		sv(2,j) = sv(2,j) - bumpdist*ryij
           		sv(3,j) = sv(3,j) - bumpdist*rzij
!       		particles are moving away, off the square shoulder potential
        	else if (coltype(i)==5) then
			wellsq=shlddia_sq(identity(i),identity(j))
           		epsave=-epsilon(1)
           		del_pe=4.0*wellsq*epsave/rmass
	   		bumpdist=smdist*dsqrt(wellsq)
	   		ratio=rmass*(-dsqrt(-del_pe+bij*bij)+bij)/(2.d0*wellsq)
           		coltype(i)=24
           		sv(1,i) = sv(1,i) + bumpdist*rxij
           		sv(2,i) = sv(2,i) + bumpdist*ryij
           		sv(3,i) = sv(3,i) + bumpdist*rzij
           		sv(1,j) = sv(1,j) - bumpdist*rxij
           		sv(2,j) = sv(2,j) - bumpdist*ryij
           		sv(3,j) = sv(3,j) - bumpdist*rzij
		else if (coltype(i)==6) then
            		wellsq=shlddia_sq(identity(i),identity(j))
           		epsave=-epsilon(1)  
           		del_pe = 4.0*wellsq*epsave/rmass
	   		bumpdist=smdist*dsqrt(wellsq)
           		if (bij*bij .gt. -del_pe) then
!    	      			particles are moving fast enough to jump into the square shoulder
	      			ratio=rmass*(dsqrt(del_pe+bij*bij)+bij)/(2.d0*wellsq)
	      			coltype(i)=25
              		sv(1,i) = sv(1,i) - bumpdist*rxij
              		sv(2,i) = sv(2,i) - bumpdist*ryij
              		sv(3,i) = sv(3,i) - bumpdist*rzij
              		sv(1,j) = sv(1,j) + bumpdist*rxij
              		sv(2,j) = sv(2,j) + bumpdist*ryij
              		sv(3,j) = sv(3,j) + bumpdist*rzij
           		else
!  	      			particles are not moving fast enough and bounce 
  	      			ratio=rmass*bij/wellsq
              		coltype(i)=26
              		sv(1,i) = sv(1,i) + bumpdist*rxij
              		sv(2,i) = sv(2,i) + bumpdist*ryij
              		sv(3,i) = sv(3,i) + bumpdist*rzij
              		sv(1,j) = sv(1,j) - bumpdist*rxij
              		sv(2,j) = sv(2,j) - bumpdist*ryij
              		sv(3,j) = sv(3,j) - bumpdist*rzij
           		endif
!       		core collision at square-shoulder diameter due to small hydrogen bond
		else if (coltype(i)==13) then
	   		wellsq=shlddia_sq(identity(i),identity(j))
	   		ratio=rmass*bij/wellsq
           		coltype(i)=27
           		bumpdist=smdist*dsqrt(wellsq)
           		sv(1,i) = sv(1,i) + bumpdist*rxij
           		sv(2,i) = sv(2,i) + bumpdist*ryij
           		sv(3,i) = sv(3,i) + bumpdist*rzij
           		sv(1,j) = sv(1,j) - bumpdist*rxij
           		sv(2,j) = sv(2,j) - bumpdist*ryij
           		sv(3,j) = sv(3,j) - bumpdist*rzij
		endif
		!LR: Changed a non-restricted else statement to a constrained one - right now i can only
		!be for species 2
        elseif(i .le. nop1+nop2) then
	 	if (coltype(i).eq.2) then
!         		bond collision
	   		blmin=ev_param(2,evcode)
        	   	ii=i-nop1-((chnnum(i)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
	   		if (evcode.eq.10) then
              		blmin=(1.d0-del_bdln(ii-3*chnln1+(chnln1*4-numbeads1)))*bdln(ii-3*chnln1+(chnln1*4-numbeads1))
	   		elseif (evcode.eq.11) then
              		blmin=(1.d0-del_blrn(ii-chnln2-3*chnln1+(chnln1*4-numbeads1)))*bl_rn(ii-chnln2-3*chnln1+(chnln1*4-numbeads1))
	   		elseif (evcode.eq.12) then
              		blmin=(1.d0-del_blrc(ii-2*chnln2-3*chnln1+(chnln1*4-numbeads1)))*bl_rc(ii-2*chnln2-3*chnln1+(chnln1*4-numbeads1))
	   		endif
	   		ratio=rmass*bij/(blmin*blmin)	   
		elseif (coltype(i).eq.3) then
!       		bond extension
	   		blmax=ev_param(3,evcode)
        	   	ii=i-nop1-((chnnum(i)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
			if (evcode.eq.10) then
              		blmax=(1.d0+del_bdln(ii-3*chnln1+(chnln1*4-numbeads1)))*bdln(ii-3*chnln1+(chnln1*4-numbeads1))
	   		elseif (evcode.eq.11) then
              		blmax=(1.d0+del_blrn(ii-chnln2-3*chnln1+(chnln1*4-numbeads1)))*bl_rn(ii-chnln2-3*chnln1+(chnln1*4-numbeads1))
	   		elseif (evcode.eq.12) then
              		blmax=(1.d0+del_blrc(ii-2*chnln2-3*chnln1+(chnln1*4-numbeads1)))*bl_rc(ii-2*chnln2-3*chnln1+(chnln1*4-numbeads1))
	   		endif
	   		ratio=rmass*bij/(blmax*blmax)
		elseif (coltype(i).eq.1) then
!       		core collision	
           		if (evcode .eq. 15) then
              		if (bptnr(i) .eq. j) then
                 			sigsq = sigma_sq(identity(i),identity(j))*ev_param(1,evcode)**2
              		else
                 			sigsq = sigma_sq(identity(i),identity(j))
              		end if
           		else
              		sigsq = sigma_sq(identity(i),identity(j))*ev_param(1,evcode)**2
              		if(evcode.ge.22.and.evcode.le.26) then
              			k=max0(identity(i),identity(j))
              			sigsq=sigsq*sqz610(evcode-21,k)**2
              		endif
           		end if
	   		ratio=rmass*bij/sigsq
        	else if (coltype(i).eq.4) then
	   		wellsq=welldia_sq(identity(i),identity(j))
           		epsave=ep_sqrt(identity(i),identity(j))
           		del_pe = 4.0*wellsq*epsave/rmass
		       if (bij*bij+del_pe.gt.0.d0) then
!          			capture
	   			ratio=rmass*(dsqrt((4.d0*wellsq*epsave/rmass)+bij*bij)+bij)/(2.d0*wellsq)
	   			coltype(i) = 20
	   			bumpdist=smdist*dsqrt(wellsq)
           			sv(1,i) = sv(1,i) - bumpdist*rxij
           			sv(2,i) = sv(2,i) - bumpdist*ryij
           			sv(3,i) = sv(3,i) - bumpdist*rzij
           			sv(1,j) = sv(1,j) + bumpdist*rxij
           			sv(2,j) = sv(2,j) + bumpdist*ryij
           			sv(3,j) = sv(3,j) + bumpdist*rzij
           		else    
! 				bounce
	      			ratio=rmass*bij/wellsq
              		coltype(i)=22
	   			bumpdist=smdist*dsqrt(wellsq)
           			sv(1,i) = sv(1,i) + bumpdist*rxij
           			sv(2,i) = sv(2,i) + bumpdist*ryij
           			sv(3,i) = sv(3,i) + bumpdist*rzij
           			sv(1,j) = sv(1,j) - bumpdist*rxij
           			sv(2,j) = sv(2,j) - bumpdist*ryij
           			sv(3,j) = sv(3,j) - bumpdist*rzij
           		endif
		else if (coltype(i)==8) then
                  	wellsq=welldia_sq(identity(i),identity(j))
           		epsave=ep_sqrt(identity(i),identity(j))
           		del_pe = 4.0*wellsq*epsave/rmass
	   		bumpdist=smdist*dsqrt(wellsq)
           		if (bij*bij.gt.del_pe) then
!          			k.e. high enough for dissociation event
	      			ratio=rmass*(-dsqrt(-del_pe+bij*bij)+bij)/(2.d0*wellsq)
              		coltype(i)=21		       
              		sv(1,i) = sv(1,i) + bumpdist*rxij
              		sv(2,i) = sv(2,i) + bumpdist*ryij
              		sv(3,i) = sv(3,i) + bumpdist*rzij
              		sv(1,j) = sv(1,j) - bumpdist*rxij
              		sv(2,j) = sv(2,j) - bumpdist*ryij  
              		sv(3,j) = sv(3,j) - bumpdist*rzij
!          			k.e. not high enough, therefore they bounce
           		else
	      			ratio=rmass*bij/wellsq
              		coltype(i)=22
              		sv(1,i) = sv(1,i) - bumpdist*rxij
              		sv(2,i) = sv(2,i) - bumpdist*ryij
              		sv(3,i) = sv(3,i) - bumpdist*rzij
              		sv(1,j) = sv(1,j) + bumpdist*rxij
              		sv(2,j) = sv(2,j) + bumpdist*ryij
              		sv(3,j) = sv(3,j) + bumpdist*rzij
           		endif
!       		core collision at square-well diameter due to overcrowding 
		else if (coltype(i)==9) then 
	   		wellsq=welldia_sq(identity(i),identity(j))
	   		ratio=rmass*bij/wellsq
           		coltype(i)=23
	   		bumpdist=smdist*dsqrt(wellsq)
           		sv(1,i) = sv(1,i) + bumpdist*rxij
           		sv(2,i) = sv(2,i) + bumpdist*ryij
           		sv(3,i) = sv(3,i) + bumpdist*rzij
           		sv(1,j) = sv(1,j) - bumpdist*rxij
           		sv(2,j) = sv(2,j) - bumpdist*ryij
           		sv(3,j) = sv(3,j) - bumpdist*rzij
!       		particles are moving away, off the square shoulder potential
        	else if (coltype(i)==5) then
			wellsq=shlddia_sq(identity(i),identity(j))
           		epsave=-epsilon(1)
           		del_pe=4.0*wellsq*epsave/rmass
	   		bumpdist=smdist*dsqrt(wellsq)
	   		ratio=rmass*(-dsqrt(-del_pe+bij*bij)+bij)/(2.d0*wellsq)
           		coltype(i)=24
           		sv(1,i) = sv(1,i) + bumpdist*rxij
           		sv(2,i) = sv(2,i) + bumpdist*ryij
           		sv(3,i) = sv(3,i) + bumpdist*rzij
           		sv(1,j) = sv(1,j) - bumpdist*rxij
           		sv(2,j) = sv(2,j) - bumpdist*ryij
           		sv(3,j) = sv(3,j) - bumpdist*rzij
		else if (coltype(i)==6) then
            		wellsq=shlddia_sq(identity(i),identity(j))
           		epsave=-epsilon(1)  
           		del_pe = 4.0*wellsq*epsave/rmass
	   		bumpdist=smdist*dsqrt(wellsq)
           		if (bij*bij .gt. -del_pe) then
!    	      			particles are moving fast enough to jump into the square shoulder
	      			ratio=rmass*(dsqrt(del_pe+bij*bij)+bij)/(2.d0*wellsq)
	      			coltype(i)=25
              		sv(1,i) = sv(1,i) - bumpdist*rxij
              		sv(2,i) = sv(2,i) - bumpdist*ryij
              		sv(3,i) = sv(3,i) - bumpdist*rzij
              		sv(1,j) = sv(1,j) + bumpdist*rxij
              		sv(2,j) = sv(2,j) + bumpdist*ryij
              		sv(3,j) = sv(3,j) + bumpdist*rzij
           		else
!  	      			particles are not moving fast enough and bounce 
  	      			ratio=rmass*bij/wellsq
              		coltype(i)=26
              		sv(1,i) = sv(1,i) + bumpdist*rxij
              		sv(2,i) = sv(2,i) + bumpdist*ryij
              		sv(3,i) = sv(3,i) + bumpdist*rzij
              		sv(1,j) = sv(1,j) - bumpdist*rxij
              		sv(2,j) = sv(2,j) - bumpdist*ryij
              		sv(3,j) = sv(3,j) - bumpdist*rzij
           		endif
!       		core collision at square-shoulder diameter due to small hydrogen bond
		else if (coltype(i)==13) then
	   		wellsq=shlddia_sq(identity(i),identity(j))
	   		ratio=rmass*bij/wellsq
           		coltype(i)=27
           		bumpdist=smdist*dsqrt(wellsq)
           		sv(1,i) = sv(1,i) + bumpdist*rxij
           		sv(2,i) = sv(2,i) + bumpdist*ryij
           		sv(3,i) = sv(3,i) + bumpdist*rzij
           		sv(1,j) = sv(1,j) - bumpdist*rxij
           		sv(2,j) = sv(2,j) - bumpdist*ryij
           		sv(3,j) = sv(3,j) - bumpdist*rzij
		endif
		!LR: This comment applies to all code between LRStart and LREnd - 
		!LR: This block is an exact duplicate of the 1st and 2nd species event collision
		! dynamics.  I know that the 3rd species cannot have all collision types present
		! but I wanted the capability there for future work.  Again, I exactly copied the 
		! 2nd species block (which was added by Dave Latshaw), I just changed the variable bounds
		! where necessary.
		!LRStart - start of 3rd species block
		
	endif

	delvx=ratio*rxij
	delvy=ratio*ryij
	delvz=ratio*rzij
	sv(4,i)=sv(4,i)-delvx/bm(i)
	sv(4,j)=sv(4,j)+delvx/bm(j)
	sv(5,i)=sv(5,i)-delvy/bm(i)
	sv(5,j)=sv(5,j)+delvy/bm(j)
	sv(6,i)=sv(6,i)-delvz/bm(i)
	sv(6,j)=sv(6,j)+delvz/bm(j)
!	rewind i,j to false positions
	sv(1,i)=sv(1,i)+delvx*tfalse/bm(i)
	sv(2,i)=sv(2,i)+delvy*tfalse/bm(i)
	sv(3,i)=sv(3,i)+delvz*tfalse/bm(i)        
	sv(1,j)=sv(1,j)-delvx*tfalse/bm(j)
	sv(2,j)=sv(2,j)-delvy*tfalse/bm(j)
	sv(3,j)=sv(3,j)-delvz*tfalse/bm(j)

	return

	end
