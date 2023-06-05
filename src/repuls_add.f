	subroutine repuls_add(i,j)

!	add nbors of hydrogen-bonded pair to the extrap_repulse list to make them
!   	appear to each other larger than they actually are

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

	ev_code(ncaj,i)=xrepuls2
    	ev_code(i,ncaj)=xrepuls1
	ev_code(i,nnjp1)=xrepuls1
     	ev_code(nnjp1,i)=xrepuls2
 	ev_code(ncai,j)=xrepuls2
    	ev_code(j,ncai)=xrepuls1
	ev_code(j,ncim1)=xrepuls1
   	ev_code(ncim1,j)=xrepuls2
	extra_repuls(i,1) = ncaj
	extra_repuls(i,2) = nnjp1
	extra_repuls(ncaj,3)=i
	extra_repuls(nnjp1,3)=i
	extra_repuls(j,1) = ncai
	extra_repuls(j,2) = ncim1
	extra_repuls(ncai,3)=j
	extra_repuls(ncim1,3)=j
    	extra_repuls(i,4)=j
    	extra_repuls(j,4)=i

 	return

      	end
