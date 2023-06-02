
	subroutine avscon

#include "def.h"

	use global
	use inputreadin

	implicit none

	real*4 xtmp(nop),ytmp(nop),ztmp(nop)
	real*8 sv_temp(6,nop)
	integer k,m,l,nframe
	character*64 filename
	save nframe
	data nframe /0/

	nframe = nframe+1

	do k = 1,nop
		m = (chnnum(k)-1)*numbeads
           l= m + (position(k-m))
           sv_temp(1,k)=sv(1,k)+sv(4,k)*tfalse
           sv_temp(2,k)=sv(2,k)+sv(5,k)*tfalse
           sv_temp(3,k)=sv(3,k)+sv(6,k)*tfalse

           xtmp(l) = sngl(sv_temp(1,k))
           ytmp(l) = sngl(sv_temp(2,k))
           ztmp(l) = sngl(sv_temp(3,k))
        enddo

        write(799) nframe, xtmp, ytmp, ztmp

	do k = 1,nop
	   sv_temp(1,k)=sv_temp(1,k)*boxl_orig
           sv_temp(2,k)=sv_temp(2,k)*boxl_orig
           sv_temp(3,k)=sv_temp(3,k)*boxl_orig
           sv_temp(4,k)=sv(4,k)
           sv_temp(5,k)=sv(5,k)
           sv_temp(6,k)=sv(6,k)
	enddo

        filename = 'results/run'//fname_digits//'.lastconfig'
        open(unit=707, file=filename)
        write(707,*) nop
        write(707,*) boxl_orig
        write(707,*) sv_temp
        close(unit=707)

        return

        end


