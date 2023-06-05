!	subroutine bond finds bond event time and type for i,j pair
!	2=bond collision, 3=bond stretch		  
	
	subroutine bond(i,j,evcode,tij,type)

#include "def.h"

	use global
	use inputreadin
	
	implicit none

	integer evcode,type
	integer i,j,ii
	real*8 tij,blmin,blmax,vxij,vyij,vzij,rxij,ryij,rzij,bij,rijsq,vijsq
	real*8 discr1,discr2	

       if (i .le. nop1) then
		ii=i-((chnnum(i)-1)*numbeads1)
!		decide how much to let bonds fluctuate with del in header.f
		blmin=ev_param(2,evcode)
		blmax=ev_param(3,evcode)
		if (evcode.eq.10) then
			blmin=(1.d0-del_bdln(ii))*bdln(ii)
			blmax=(1.d0+del_bdln(ii))*bdln(ii)
		elseif (evcode.eq.11) then
			blmin=(1.d0-del_blrn(ii-chnln1))*bl_rn(ii-chnln1)
			blmax=(1.d0+del_blrn(ii-chnln1))*bl_rn(ii-chnln1)
		elseif (evcode.eq.12) then
			blmin=(1.d0-del_blrc(ii-2*chnln1))*bl_rc(ii-2*chnln1)
			blmax=(1.d0+del_blrc(ii-2*chnln1))*bl_rc(ii-2*chnln1)
		endif
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
!       		then bond collision and bond stretch are both possible
	   		discr1=bij*bij-vijsq*(rijsq-blmin*blmin)
	   		if (discr1.gt.0.d0) then
!             		bond collision will occur at bl-delta
	      			tij=(-bij-dsqrt(discr1))/vijsq
              		type=2
	   		else
!          			bond extension will occur at bl+delta
	       		discr2=bij*bij-vijsq*(rijsq-blmax*blmax)
	       		if (discr2.gt.0.d0) then
	          			tij=(-bij+dsqrt(discr2))/vijsq
	          			type=3
	       		endif
	   		endif
		else
!			bond extension will occur at bl+delta
			discr2=bij*bij-vijsq*(rijsq-blmax*blmax)
			if (discr2.gt.0.d0) then
				tij=-(rijsq-blmax*blmax)/(dsqrt(discr2)+bij)
				type=3
			endif
		endif
       else
       	ii=i-nop1-((chnnum(i)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
!		decide how much to let bonds fluctuate with del in header.f
		blmin=ev_param(2,evcode)
		blmax=ev_param(3,evcode)
		if (evcode.eq.10) then
           		blmin=(1.d0-del_bdln(ii-3*chnln1+(chnln1*4-numbeads1)))*bdln(ii-3*chnln1+(chnln1*4-numbeads1))
           		blmax=(1.d0+del_bdln(ii-3*chnln1+(chnln1*4-numbeads1)))*bdln(ii-3*chnln1+(chnln1*4-numbeads1))
		elseif (evcode.eq.11) then
          		blmin=(1.d0-del_blrn(ii-chnln2-3*chnln1+(chnln1*4-numbeads1)))*bl_rn(ii-chnln2-3*chnln1+(chnln1*4-numbeads1))
           		blmax=(1.d0+del_blrn(ii-chnln2-3*chnln1+(chnln1*4-numbeads1)))*bl_rn(ii-chnln2-3*chnln1+(chnln1*4-numbeads1))
		elseif (evcode.eq.12) then
           		blmin=(1.d0-del_blrc(ii-2*chnln2-3*chnln1+(chnln1*4-numbeads1)))*bl_rc(ii-2*chnln2-3*chnln1+(chnln1*4-numbeads1))
           		blmax=(1.d0+del_blrc(ii-2*chnln2-3*chnln1+(chnln1*4-numbeads1)))*bl_rc(ii-2*chnln2-3*chnln1+(chnln1*4-numbeads1))
		endif
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
!       		then bond collision and bond stretch are both possible
	   		discr1=bij*bij-vijsq*(rijsq-blmin*blmin)
	   		if (discr1.gt.0.d0) then
!             		bond collision will occur at bl-delta
	      			tij=(-bij-dsqrt(discr1))/vijsq
              		type=2
	   		else
!          			bond extension will occur at bl+delta
	       		discr2=bij*bij-vijsq*(rijsq-blmax*blmax)
	       		if (discr2.gt.0.d0) then
	          			tij=(-bij+dsqrt(discr2))/vijsq
	          			type=3
	       		end if
	   		endif
		else
!       		bond extension will occur at bl+delta
	   		discr2=bij*bij-vijsq*(rijsq-blmax*blmax)
	   		if (discr2.gt.0.d0) then
	      			tij=-(rijsq-blmax*blmax)/(dsqrt(discr2)+bij)
	      			type=3
	   		end if
		endif
       endif

	return

	end
	
