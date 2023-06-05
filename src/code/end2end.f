	subroutine end_to_end(e2e_avg)

#include "def.h"

	use global
	use inputreadin
	implicit none

      	real*8 e2ex,e2ey,e2ez,e2e,e2e_avg,vxij,vyij,vzij
	integer l,ll
      
	e2e_avg = 0.d0

	do l=1,(nop1/numbeads1)
		ll=(l-1)*numbeads1   
		vxij=sv(4,ll+chnln1+1)-sv(4,ll+3*chnln1)
           	vyij=sv(5,ll+chnln1+1)-sv(5,ll+3*chnln1)
           	vzij=sv(6,ll+chnln1+1)-sv(6,ll+3*chnln1)
           	e2ex=sv(1,ll+chnln1+1) - sv(1,ll+3*chnln1) + vxij*tfalse
          	e2ey=sv(2,ll+chnln1+1) - sv(2,ll+3*chnln1) + vyij*tfalse
           	e2ez=sv(3,ll+chnln1+1) - sv(3,ll+3*chnln1) + vzij*tfalse
	  	e2ex=e2ex-dnint(e2ex)
	   	e2ey=e2ey-dnint(e2ey)
	   	e2ez=e2ez-dnint(e2ez)
	   	e2e=dsqrt(e2ex**2+e2ey**2+e2ez**2)
	   	e2e_avg = e2e_avg + e2e
	end do

	do l=(nop1/numbeads1)+1,numchains
		ll=nop1+(l-nop1/numbeads1-1)*numbeads2   
           	vxij=sv(4,ll+chnln2+1)-sv(4,ll+3*chnln2)
           	vyij=sv(5,ll+chnln2+1)-sv(5,ll+3*chnln2)
           	vzij=sv(6,ll+chnln2+1)-sv(6,ll+3*chnln2)
           	e2ex=sv(1,ll+chnln2+1) - sv(1,ll+3*chnln2) + vxij*tfalse
           	e2ey=sv(2,ll+chnln2+1) - sv(2,ll+3*chnln2) + vyij*tfalse
           	e2ez=sv(3,ll+chnln2+1) - sv(3,ll+3*chnln2) + vzij*tfalse
	   	e2ex=e2ex-dnint(e2ex)
	   	e2ey=e2ey-dnint(e2ey)
	   	e2ez=e2ez-dnint(e2ez)
	  	e2e=dsqrt(e2ex**2+e2ey**2+e2ez**2)
	   	e2e_avg = e2e_avg + e2e
	end do
	
	e2e_avg = e2e_avg/dble(numchains)*boxl_orig

	return

	end




