	subroutine scale_down()

#include "def.h"

   	use global
	use inputreadin

      	implicit none

      	integer k,kk

!     	scale positions, bond lengths, bead sizes, and well sizes 
!     	so that boxl=1

    	write(fileout,*)' '
   	write(fileout,*)'scaling, original boxl=',boxl
   	boxl_orig=boxl

#ifndef runr
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
     	do k=1,(noptotal)
      		sv(1,k)=sv(1,k)/boxl
       	sv(2,k)=sv(2,k)/boxl
         	sv(3,k)=sv(3,k)/boxl
      	enddo
#endif

    	do k=1,28
	 	sigma(k)=sigma(k)/boxl_orig
	 	welldia(k)=welldia(k)/boxl_orig
      	enddo

      	shder_dist1 = 5.00d0/boxl_orig
      	shder_dist2 = 4.74d0/boxl_orig
      	shder_dist3 = 4.86d0/boxl_orig
      	shder_dist4 = 4.83d0/boxl_orig

      	do k=1,28
         	do kk = 1,28
	    		sigma_sq(k,kk) = 0.25d0*(sigma(k)+sigma(kk))**2
	    		sigma_2b(k,kk) = 0.50d0*(sigma(k)+sigma(kk))
	    		welldia_sq(k,kk) = 0.25d0*(welldia(k)+welldia(kk))**2
	    		ep_sqrt(k,kk) = dsqrt(epsilon(k)*epsilon(kk))
	    		shlddia_sq(k,kk) = shder_dist4**2
         	enddo
   	enddo

      	do k = 9,28
         	do kk = 9,28
            		ep_sqrt(k,kk) = epsilon(1)*ep(k,kk)
            		sigma_sq(k,kk) = bds(k,kk)**2/boxl_orig**2
	    		sigma_2b(k,kk) = bds(k,kk)/boxl_orig
            		welldia_sq(k,kk) = wel(k,kk)**2/boxl_orig**2
         	enddo
   	enddo
	    
      	shlddia_sq(1,2) = shder_dist1**2
      	shlddia_sq(2,1) = shder_dist1**2
      	shlddia_sq(5,2) = shder_dist1**2
      	shlddia_sq(2,5) = shder_dist1**2
      	shlddia_sq(1,1) = shder_dist2**2
      	shlddia_sq(5,5) = shder_dist2**2
      	shlddia_sq(1,5) = shder_dist2**2
      	shlddia_sq(5,1) = shder_dist2**2
      	shlddia_sq(2,4) = shder_dist3**2
      	shlddia_sq(4,2) = shder_dist3**2
      	shlddia_sq(2,8) = shder_dist3**2
      	shlddia_sq(8,2) = shder_dist3**2

      	do k=1,chnln1
         	bdln(k)=bdln(k)/boxl_orig
         	bl_rn(k)=bl_rn(k)/boxl_orig
         	bl_rc(k)=bl_rc(k)/boxl_orig
      	enddo

      	do k=1,chnln2
         	bdln(chnln1+k)=bdln(chnln1+k)/boxl_orig
         	bl_rn(chnln1+k)=bl_rn(chnln1+k)/boxl_orig
         	bl_rc(chnln1+k)=bl_rc(chnln1+k)/boxl_orig
      	enddo

      	boxl=1.d0
      	write(fileout,*)'done scaling, boxl=',boxl
      	half=boxl/2.d0

      	return

      	end
