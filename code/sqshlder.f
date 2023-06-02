!	subroutine sqwel finds square shoulder event time and type for i,j
!      pair; 4=sw capture, 5=sw dissociation, 6=sw bounce
!      only eligible side chain pairs make it into this subroutine

	subroutine sqshlder(i,j,evcode,tij,type)

#include "def.h"

       use global
	use inputreadin
	
	implicit none

	integer evcode,i,j,type
       real*8 vxij,vyij,vzij,rxij,ryij,rzij,bij,rijsq,vijsq
       real*8 diff,corediscr,shlddiscr,b_new,tij

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
	rijsq=rxij*rxij+ryij*ryij+rzij*rzij
	vijsq=vxij*vxij+vyij*vyij+vzij*vzij
	if (bij.lt.0.d0) then
!      	then i and j are moving together
	    	diff=rijsq-shlddia_sq(identity(i),identity(j))
	    	if (diff.lt.0.d0) then
!              	then shoulders are already overlapping
	        	corediscr=bij*bij-vijsq*(rijsq-sigma_sq(identity(i),identity(j)))
	        	if (corediscr.gt.0.d0) then
!                  		then a core collision will occur
		     		tij=(-bij-dsqrt(corediscr))/vijsq
		     		type=1
               	else 
                   		shlddiscr=bij*bij-vijsq*diff
		     		tij=(-bij+dsqrt(shlddiscr))/vijsq
		     		type=10
	        	endif
	    	else
!              	the shoulders are not yet overlapping
               	shlddiscr=bij*bij-vijsq*diff
	        	if (shlddiscr.gt.0.d0) then
!                  		then a capture will occur
		     		tij=(-bij-dsqrt(shlddiscr))/vijsq
		     		type=12
	        	endif
	    	endif
	else
!          	i and j are moving away from each other
           	diff=rijsq-shlddia_sq(identity(i),identity(j))
	    	if (diff.lt.0.d0) then
!              	then wells are already overlapping and will reach discontinuity at tij
               	shlddiscr=bij*bij-vijsq*diff
	        	tij=(-bij+dsqrt(shlddiscr))/vijsq
	        	type=10
	    	endif
	endif

	return

	end
	
