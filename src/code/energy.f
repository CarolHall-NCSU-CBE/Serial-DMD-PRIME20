      	subroutine energy(ered,tred,sumvel,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij)
    
#include "def.h"

      	use global
	use inputreadin

      	implicit none

      	real*8 sumeps,eps_hb,eps_hh,welli,epsi,wellj,epsj,wellsq,vxij,vyij,vzij
      	real*8 rxij,ryij,rzij,rijsq,ered,tred,sumvel
     	integer i,j,k,sum_hb,hb_ii,hb_ij,sum_hh,hh_ii,hh_ij,hb_alpha
	real*8 ehh_ii,ehh_ij,ep_depth,sum_ehh  !! by mookyung

      	hb_alpha = 0
      	hb_ii = 0
      	hb_ij = 0
      	hh_ii = 0
      	hh_ij = 0
      	ehh_ii = 0.0d0
      	ehh_ij = 0.0d0
!     	only want to consider beads that are "allowed" to have square
!	well interactions as per events.f and eventredo.f

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do i=1,(noptotal)-1
	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
		do j=i+1,(noptotal)
!           		count the number of hb's using the array bptnr
	    		if (bptnr(i) == j) then
	       		if (chnnum(i) == chnnum(j)) then
	          			hb_ii = hb_ii + 1
	       		else
	          			hb_ij = hb_ij + 1
	       		end if
!           			count the number of hydrophobic side-chain overlaps
            		else if (ev_code(j,i).eq.16) then
               		vxij=sv(4,i)-sv(4,j)
               		vyij=sv(5,i)-sv(5,j)	
               		vzij=sv(6,i)-sv(6,j)
               		rxij=sv(1,i)-sv(1,j)+vxij*tfalse
               		ryij=sv(2,i)-sv(2,j)+vyij*tfalse
               		rzij=sv(3,i)-sv(3,j)+vzij*tfalse
               		rxij=rxij-dnint(rxij)
               		ryij=ryij-dnint(ryij)
               		rzij=rzij-dnint(rzij)
               		rijsq=rxij*rxij+ryij*ryij+rzij*rzij
               		wellsq=welldia_sq(identity(i),identity(j))
               		ep_depth=ep_sqrt(identity(i),identity(j))
               		if (rijsq.le.wellsq) then		  
	          			if (chnnum(i) == chnnum(j)) then
                     			ehh_ii = ehh_ii + ep_depth
                 			else 
                     			ehh_ij = ehh_ij + ep_depth
                  			end if
               		endif
            		endif
         	enddo
      enddo     

	do k=1,(nop1/numbeads1)
		do i=(k-1)*numbeads1+chnln1+5,(k-1)*numbeads1+2*chnln1
			j=i+chnln1-4
			if (bptnr(i) == j) hb_alpha = hb_alpha + 1
		enddo
	enddo

	do k=1,(nop2/numbeads2)
		do i=nop1+(k-1)*numbeads2+chnln2+5,nop1+(k-1)*numbeads2+2*chnln2
			j=i+chnln2-4
			if (bptnr(i) == j) hb_alpha = hb_alpha + 1
		enddo
	enddo

	eps_hb=ep_sqrt(5,8)
      	eps_hh=ep_sqrt(20,20)      
      	sum_hb = hb_ii + hb_ij
      	sum_ehh = ehh_ii + ehh_ij
      	sumeps = -(sum_hb*eps_hb+sum_ehh)
!     	sum velocities^2
      	sumvel=0.d0

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do i=1,(noptotal)
#ifdef glycine 
!		to not include glycine side-chains
		if ( i .le. (nop1/numeabds1) then
			if (bdln(i-((chnnum(i)-1)*numbeads1)).lt.99) then
		else
			if (bdln(i-nop1-((chnnum(i)-(nop1/numbeads1)-1)*numbeads2)+numbeads1)).lt.99) then
#endif
		sumvel=sumvel+bm(i)*(sv(4,i)*sv(4,i)+sv(5,i)*sv(5,i)+sv(6,i)*sv(6,i))
#ifdef glycine
		endif
#endif
	enddo

!	calculate reduced energy (ered=total energy )
      	ered=0.5d0*sumvel+sumeps
!     	calculate reduced temperature (tred=k*temperature)
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	tred=sumvel/3.d0/dble((noptotal))
  
      	return

      	end
