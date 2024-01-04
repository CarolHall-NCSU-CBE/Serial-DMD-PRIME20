	subroutine scale_up()

#include "def.h"

      	use global
	use inputreadin

      	implicit none

      	integer k

!     	scale positions, bead sizes, and well sizes back to original
      	write(6,*)' '
      	write(6,*)'scaling, boxl was=',boxl
      	boxl=boxl_orig

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k=1,(noptotal)
         	sv(1,k)=sv(1,k)*boxl_orig
         	sv(2,k)=sv(2,k)*boxl_orig
         	sv(3,k)=sv(3,k)*boxl_orig
      	enddo

      	do k=1,28
		sigma(k)=sigma(k)*boxl_orig
		welldia(k)=welldia(k)*boxl_orig
      	enddo
	
      	do k=1,chnln1
         	bdln(k)=bdln(k)*boxl_orig
         	bl_rn(k)=bl_rn(k)*boxl_orig
         	bl_rc(k)=bl_rc(k)*boxl_orig
      	enddo

      	do k=1,chnln2
         	bdln(chnln1+k)=bdln(chnln1+k)*boxl_orig
         	bl_rn(chnln1+k)=bl_rn(chnln1+k)*boxl_orig
         	bl_rc(chnln1+k)=bl_rc(chnln1+k)*boxl_orig
      	enddo

      	write(6,*)'done scaling, boxl=',boxl
         
      	return

      	end
