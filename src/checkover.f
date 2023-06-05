!	subroutine checkover checks entire configuration for particle 
!	overlaps 
!	only consider non-bonded pairs where neither i nor j is a gly

	subroutine checkover(over)
 
#include "def.h"

	use global
	use inputreadin

	implicit none

	logical over
	integer i,j,k,ii,jj,evcode
	real*8 sigsq,rxij,ryij,rzij,rijsq,vxij,vyij,vzij,blmin,blmax
	real*8 rijsq_min,rijsq_max

	over=.false.

	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do i=1,(noptotal)-1
		if (i .le. nop1) then
	   		ii=i-(chnnum(i)-1)*numbeads1
	   		!LR: Changed an open-ended else to a constrained one, since i could now be greater than nop1+nop2
          	elseif (i .le. nop1+nop2) then
          		ii=i-nop1-((chnnum(i)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
          	endif
          	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	      	do j=i+1,(noptotal)
			if (j .le. nop1) then
		 		jj=j-(chnnum(j)-1)*numbeads1
		 		!LR: Changed an open-ended else to a constrained one, since j could now be greater than nop1+nop2
               	elseif (j .le. nop1+nop2) then
               		jj=j-nop1-((chnnum(j)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
               	endif
	         	evcode = ev_code(j,i)
		    	if ((evcode.eq.1).or.(evcode.ge.15)) then
	              	sigsq=sigma_sq(identity(i),identity(j))*ev_param(1,evcode)**2
                       	if(evcode.ge.22.and.evcode.le.26) then
                       		k=max0(identity(i),identity(j))
                       		sigsq=sigsq*sqz610(evcode-21,k)**2
                       	endif
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
	               	rijsq=rijsq*1.0000000001d0
	               	if (rijsq.le.sigsq) then		
	                  		over=.true.
	                  		write(fileout,*)'particles ',i,' and ',j,' overlap'
	                  		write(fileout,*)'ev_code(i,j)=',evcode
	                  		write(fileout,*)'identity(i)',identity(i)
	                  		write(fileout,*)'identity(j)',identity(j)
	                  		write(fileout,*)'factor=',ev_param(1,evcode)
	                  		write(fileout,*)'rijsq=',rijsq*boxl_orig*boxl_orig
	                  		write(fileout,*)'rij=',dsqrt(rijsq)*boxl_orig
	                  		write(fileout,*)'sigsq=',sigsq*boxl_orig*boxl_orig
	                  		write(fileout,*)'sig=',dsqrt(sigsq)*boxl_orig
	               	endif
		    	elseif ((evcode.ge.4).and.(evcode.le.12)) then
				if (chaptype .eq. 1) then
		      			blmin=ev_param(2,evcode)
		       		blmax=ev_param(3,evcode)
                       		if (i .le. nop1) then
                       			if (evcode.eq.10) then
                          				blmin=(1.d0-del_bdln(ii))*bdln(ii)
                          				blmax=(1.d0+del_bdln(ii))*bdln(ii)
                       			elseif (evcode.eq.11) then
                          				blmin=(1.d0-del_blrn(ii-chnln1))*bl_rn(ii-chnln1)
                          				blmax=(1.d0+del_blrn(ii-chnln1))*bl_rn(ii-chnln1)
                       			elseif (evcode.eq.12) then
                         				blmin=(1.d0-del_blrc(ii-2*chnln1))*bl_rc(ii-2*chnln1)
                          				blmax=(1.d0+del_blrc(ii-2*chnln1))*bl_rc(ii-2*chnln1)
                       			endif
                       		elseif (i .le. (nop1+nop2)) then
                       			if (evcode.eq.10) then
                          				blmin=(1.d0-del_bdln(ii-3*chnln1+(chnln1*4-numbeads1)))*bdln(ii-3*chnln1+(chnln1*4-numbeads1))
                          				blmax=(1.d0+del_bdln(ii-3*chnln1+(chnln1*4-numbeads1)))*bdln(ii-3*chnln1+(chnln1*4-numbeads1))
                       			elseif (evcode.eq.11) then
                          				blmin=(1.d0-del_blrn(ii-chnln2-3*chnln1+(chnln1*4-numbeads1)))*bl_rn(ii-chnln2-3*chnln1+(chnln1*4-numbeads1))
                          				blmax=(1.d0+del_blrn(ii-chnln2-3*chnln1+(chnln1*4-numbeads1)))*bl_rn(ii-chnln2-3*chnln1+(chnln1*4-numbeads1))
                       			elseif (evcode.eq.12) then
                          				blmin=(1.d0-del_blrc(ii-2*chnln2-3*chnln1+(chnln1*4-numbeads1)))*bl_rc(ii-2*chnln2-3*chnln1+(chnln1*4-numbeads1))
                          				blmax=(1.d0+del_blrc(ii-2*chnln2-3*chnln1+(chnln1*4-numbeads1)))*bl_rc(ii-2*chnln2-3*chnln1+(chnln1*4-numbeads1))
                       			endif
					endif

		       		blmin=blmin**2
	              		blmax=blmax**2
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
		       		rijsq_min=rijsq*1.0000000001d0	
		       		rijsq_max=rijsq*0.9999999999d0
	
		       		if (rijsq_max .gt. blmax) then 
						over=.true.
                          			write(fileout,*)'bond between particles ',i,' and ',j,' are too long'
                          			write(fileout,*)'ev_code(i,j)=',evcode
                          			write(fileout,*)'identity(i)',identity(i)
                          			write(fileout,*)'identity(j)',identity(j)
                          			write(fileout,*)'bond=',dsqrt(blmax)*boxl_orig
                          			write(fileout,*)'rij=',dsqrt(rijsq)*boxl_orig
		       		end if

					if (rijsq_min .lt. blmin) then
                          			over=.true.
                          			write(fileout,*)'bond between particles ',i,' and ',j,' are too short'
                          			write(fileout,*)'ev_code(i,j)=',evcode
                          			write(fileout,*)'identity(i)',identity(i)
                         			write(fileout,*)'identity(j)',identity(j)
                          			write(fileout,*)'bond=',dsqrt(blmin)*boxl_orig
                          			write(fileout,*)'rij=',dsqrt(rijsq)*boxl_orig
					end if
				endif
			endif	
		enddo
	enddo

	return

	end

