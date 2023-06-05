	subroutine sum(hb,hh)

#include "def.h"

      	use global
	use inputreadin

      	implicit none

      	integer hb,hh
        
      	if (coll .le. 2*quarter) then
		hb_sum(2) = hb_sum(2) + hb
		hh_sum(2) = hh_sum(2) + hh
		qt_count(2) = qt_count(2) + 1
	elseif (coll .le. 3*quarter) then
		hb_sum(3) = hb_sum(3) + hb
		hh_sum(3) = hh_sum(3) + hh
		qt_count(3) = qt_count(3) + 1
	else
		hb_sum(4) = hb_sum(4) + hb
		hh_sum(4) = hh_sum(4) + hh
		qt_count(4) = qt_count(4) + 1
	endif

      	return
      
      	end

!*********************************************************************************

	subroutine average(per_hb,per_hh)

#include "def.h"

	use global

      	implicit none

      	integer k
      	real per_hb,per_hh
      	real av_hb(4),av_hh(4),sqsum_hb,sumsq_hb,sqsum_hh,sumsq_hh
      	character*64 filename
      
      	sqsum_hb = 0.0
      	sumsq_hb = 0.0
      	sqsum_hh = 0.0
      	sumsq_hh = 0.0

      	do k = 2, 4
		if (qt_count(k) .ne. 0) then
			av_hb(k) = real(hb_sum(k))/qt_count(k)
			av_hh(k) = real(hh_sum(k))/qt_count(k)
	 	else
	    		av_hb(k) = 0.0
	    		av_hh(k) = 0.0
	 	end if
         	sqsum_hb = sqsum_hb + av_hb(k)**2
         	sqsum_hh = sqsum_hh + av_hh(k)**2
	 	sumsq_hb = sumsq_hb + av_hb(k)
	 	sumsq_hh = sumsq_hh + av_hh(k)
	end do  
    
      	sumsq_hb = sumsq_hb**2
      	sumsq_hh = sumsq_hh**2
      
	if (sumsq_hb .ne. 0.0) then
		per_hb = sqrt(abs(3.0*sqsum_hb - sumsq_hb)/(2.0/3.0*sumsq_hb))*100.0
	else
	 	per_hb = 0.0
      	endif

      	if (sumsq_hh .ne. 0.0) then
         	per_hh = sqrt(abs(3.0*sqsum_hh - sumsq_hh)/(2.0/3.0*sumsq_hh))*100.0
      	else
	 	per_hh = 0.0
      	endif
	
      	write(fileout,*) 'av_hb', av_hb(2), av_hb(3), av_hb(4), per_hb
      	write(fileout,*) 'av_hh', av_hh(2), av_hh(3), av_hh(4), per_hh

      	if (av_hb(4) .le. 10.0) per_hb = 0.0
      	if (av_hh(4) .le. 10.0) per_hh = 0.0

      	return

      	end

!*********************************************************************************

      	subroutine doubling

#include "def.h"

      	use global

      	implicit none

      	hb_sum(2) = hb_sum(3) + hb_sum(4)
      	hh_sum(2) = hh_sum(3) + hh_sum(4)
      	qt_count(2) = qt_count(3) + qt_count(4)
      	hb_sum(3) = 0
      	hb_sum(4) = 0
      	hh_sum(3) = 0
      	hh_sum(4) = 0
      	qt_count(3) = 0
      	qt_count(4) = 0

      	return

      	end

!*********************************************************************************
