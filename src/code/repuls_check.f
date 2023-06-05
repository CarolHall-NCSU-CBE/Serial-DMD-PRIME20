	subroutine repuls_check(i,j,rating)

!	add nbors of hydrogen-bonded pair to the extrap_repulse list to make them
!   	appear to each other larger than they actually are

#include "def.h"

    	use global
	use inputreadin

	implicit none

  	real*8 vxij,vyij,vzij,rxij,ryij,rzij,rating
    	real*8 d1,d2,d3,d4
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

     	vxij=sv(4,i)-sv(4,ncaj)
     	vyij=sv(5,i)-sv(5,ncaj)
   	vzij=sv(6,i)-sv(6,ncaj) 
    	rxij=sv(1,i)-sv(1,ncaj)+vxij*tfalse
	ryij=sv(2,i)-sv(2,ncaj)+vyij*tfalse
	rzij=sv(3,i)-sv(3,ncaj)+vzij*tfalse                      
	rxij=rxij-dnint(rxij)
	ryij=ryij-dnint(ryij)
	rzij=rzij-dnint(rzij)
	d1=(rxij)**2+(ryij)**2+(rzij)**2
	d1=dsqrt(d1)
       vxij=sv(4,i)-sv(4,nnjp1)
	vyij=sv(5,i)-sv(5,nnjp1)
	vzij=sv(6,i)-sv(6,nnjp1)
	rxij=sv(1,i)-sv(1,nnjp1)+vxij*tfalse
	ryij=sv(2,i)-sv(2,nnjp1)+vyij*tfalse
	rzij=sv(3,i)-sv(3,nnjp1)+vzij*tfalse
	rxij=rxij-dnint(rxij)
	ryij=ryij-dnint(ryij)
	rzij=rzij-dnint(rzij)
	d2=(rxij)**2+(ryij)**2+(rzij)**2
	d2=dsqrt(d2)
	vxij=sv(4,j)-sv(4,ncai)
	vyij=sv(5,j)-sv(5,ncai)
	vzij=sv(6,j)-sv(6,ncai)
	rxij=sv(1,j)-sv(1,ncai)+vxij*tfalse
	ryij=sv(2,j)-sv(2,ncai)+vyij*tfalse
	rzij=sv(3,j)-sv(3,ncai)+vzij*tfalse
	rxij=rxij-dnint(rxij)
	ryij=ryij-dnint(ryij)
	rzij=rzij-dnint(rzij)
	d3=(rxij)**2+(ryij)**2+(rzij)**2
	d3=dsqrt(d3)
	vxij=sv(4,j)-sv(4,ncim1)
	vyij=sv(5,j)-sv(5,ncim1)
	vzij=sv(6,j)-sv(6,ncim1)
	rxij=sv(1,j)-sv(1,ncim1)+vxij*tfalse
	ryij=sv(2,j)-sv(2,ncim1)+vyij*tfalse
	rzij=sv(3,j)-sv(3,ncim1)+vzij*tfalse
	rxij=rxij-dnint(rxij)
	ryij=ryij-dnint(ryij)
	rzij=rzij-dnint(rzij)
       d4=(rxij)**2+(ryij)**2+(rzij)**2
       d4=dsqrt(d4)

	if ((d1 > shder_dist1) .and. (d2 > shder_dist2) .and. (d3 > shder_dist3) .and. (d4 > shder_dist4)) then
   		rating = 10.d0
     	else
  		rating = 11.d0
    	endif

   	return

    	end
