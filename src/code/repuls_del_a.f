	subroutine repuls_del_a(i,j)

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

    	ev_code(ncaj,i)=1
      	ev_code(i,ncaj)=1
     	ev_code(i,nnjp1)=1
     	ev_code(nnjp1,i)=1
      	ev_code(ncai,j)=1
     	ev_code(j,ncai)=1
      	ev_code(ncim1,j)=1
      	ev_code(j,ncim1)=1

    	return

      	end


