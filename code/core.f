!	subroutine core finds core collision time for i,j pair
			  
	subroutine core (i,j,evcode,tij,type)

#include "def.h"

	use global
	use inputreadin
	
	implicit none

	integer evcode,type,i,j,k
	real*8 vxij,vyij,vzij,rxij,ryij,rzij,bij,sigsq,rijsq,vijsq,discr,tij
	
	vxij=sv(4,i)-sv(4,j)
	vyij=sv(5,i)-sv(5,j)
	vzij=sv(6,i)-sv(6,j)
	rxij=sv(1,i)-sv(1,j)+vxij*tfalse
	ryij=sv(2,i)-sv(2,j)+vyij*tfalse
    	rzij=sv(3,i)-sv(3,j)+vzij*tfalse
    	rxij=rxij-dnint(rxij)
   	ryij=ryij-dnint(ryij)
    	rzij=rzij-dnint(rzij)
	bij=rxij*vxij+ryij*vyij+rzij*vzij

	if (bij.lt.0.d0) then
!		then i and j are moving together
	   	sigsq=sigma_sq(identity(i),identity(j))*ev_param(1,evcode)**2
           	if(evcode.ge.22.and.evcode.le.26) then
           		k=max0(identity(i),identity(j))
           		sigsq=sigsq*sqz610(evcode-21,k)**2
           	endif
           	rijsq=rxij*rxij+ryij*ryij+rzij*rzij
           	vijsq=vxij*vxij+vyij*vyij+vzij*vzij
	   	discr=bij*bij-vijsq*(rijsq-sigsq)
	   	if (discr.gt.0.d0) then
!	      		then a collision occurs
	      		tij=(-bij-dsqrt(discr))/vijsq
	      		type=1
	   	endif  
	endif
	
	return

	end
	















