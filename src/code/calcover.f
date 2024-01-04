	subroutine calcover(i,j,evcode,over)

#include "def.h"

	use global
	use inputreadin
	
	implicit none

!	only i,j pairs that are not allowed to overlap should make it 
!	this far in loop
	   
	logical over
	integer i,j,k,evcode
	real*8 sigsq,rxij,ryij,rzij,rijsq

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
	rijsq=rijsq*1.000000000001d0

	if (rijsq.le.sigsq) then		
		over=.true.
		write(6,*)'particles ',i,' and ',j,' overlap'
		write(6,*)'ev_code(i,j)=',evcode
		write(6,*)identity(i)
		write(6,*)identity(j)
		write(6,*)'factor=',ev_param(1,evcode)
		write(6,*)'rijsq=',rijsq*boxl_orig*boxl_orig
		write(6,*)'rij=',dsqrt(rijsq)*boxl_orig
		write(6,*)'sigsq=',sigsq*boxl_orig*boxl_orig
		write(6,*)'sig=',dsqrt(sigsq)*boxl_orig
	endif

	return

	end

