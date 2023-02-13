	subroutine config()

#include "def.h"

	use global
	use allinone

	implicit none

    !LR: Changed a hardcoded 2-species variable reference to a noptotal variable	
	real*8 xtmp(noptotal),ytmp(noptotal),ztmp(noptotal)
	integer k
	character*64 filename

	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do k = 1,(noptotal)
		xtmp(k)=sv(1,k)+sv(4,k)*tfalse
		ytmp(k)=sv(2,k)+sv(5,k)*tfalse
		ztmp(k)=sv(3,k)+sv(6,k)*tfalse
	enddo

	filename = 'results/run'//fname_digits//'.config'
	open(unit=800,file=filename,status='unknown', form='unformatted',position='append')
	write(800) coll, t+tfalse,xtmp, ytmp, ztmp
	close(unit=800)
	filename = 'results/run'//fname_digits//'.bptnr'
	open(unit=801,file=filename,status='unknown', form='unformatted',position='append')
	write(801) coll, bptnr
	close(unit=801)

	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do k = 1,(noptotal)
		xtmp(k)=sv(4,k)
		ytmp(k)=sv(5,k)
		ztmp(k)=sv(6,k)
	enddo

	filename = 'results/run'//fname_digits//'.lastvel'
	open(unit=802,file=filename,form='unformatted')
	write(802) coll, xtmp, ytmp, ztmp
	close(unit=802)
	
        return

        end
