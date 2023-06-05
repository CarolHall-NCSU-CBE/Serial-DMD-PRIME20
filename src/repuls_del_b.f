	subroutine repuls_del_b(i,j)

!   	delete nbors of the just broken-up hydrogen-bonded pair from
!	the extrap_repulse list 

#include "def.h"

     	use global
	use inputreadin

	implicit none

	integer i,j

      	if (i .le. nop1) then
     		ncim1=i+chnln1-1
           	ncai=i-chnln1
     	else
           	ncim1=i+chnln2-1
           	ncai=i-chnln2
      	endif

     	if (j .le. nop1) then
           	ncaj=j-2*chnln1
           	nnjp1=j-chnln1+1
       else
           	ncaj=j-2*chnln2
           	nnjp1=j-chnln2+1
    	endif

      	extra_repuls(i,1) = 0 
      	extra_repuls(i,2) = 0
     	extra_repuls(ncaj,3) = 0   
     	extra_repuls(nnjp1,3) = 0
    	extra_repuls(j,1) = 0
     	extra_repuls(j,2) = 0
    	extra_repuls(ncai,3) = 0   
     	extra_repuls(ncim1,3) = 0  
    	extra_repuls(i,4)=0
     	extra_repuls(j,4)=0

    	return

     	end


