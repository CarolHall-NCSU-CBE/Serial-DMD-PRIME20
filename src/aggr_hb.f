	subroutine calcaggregates_hb(num_aggs_hb)

#include "def.h"

	use global
	use inputreadin

	implicit none

	integer i, j, n, m, p, num_aggs_hb

	do i = 1, numchains/2
		do j = 0, numchains
			newaggregate_hb(i,j) = 0
		end do
	end do

	do i = 1, numchains - 1
		do j = i+1, numchains
			if (matrix_hb(i,j) .ne. 0) then
				call addtoaggregate_hb(i,j) 
			endif
		enddo
	enddo

	num_aggs_hb = 0

	do n = 1, numchains/2
		if (newaggregate_hb(n,0) .gt. 0) then
			num_aggs_hb = num_aggs_hb + 1
			aggregate_hb(num_aggs_hb,0)=newaggregate_hb(n,0)
			do m = 1, newaggregate_hb(n,0)
				aggregate_hb(num_aggs_hb,m) = newaggregate_hb(n,m)
			end do
		endif
	end do

	return

	end 
		    
!*********************************************************************************************
		
	subroutine addtoaggregate_hb(i,j)

	use global

	implicit none

	integer i, j, n, m, indexi, indexj

	indexi = 0
	indexj = 0

	do n = 1, numchains/2
		do m = 1, newaggregate_hb(n,0)
			if (newaggregate_hb(n,m) .eq. i) then
				indexi = n
			else if (newaggregate_hb(n,m) .eq. j) then
				indexj = n
			endif
		end do
	end do

	if (indexi .eq. indexj) then
		if (indexi .eq. 0) then
			n = 1
			do while (newaggregate_hb(n,0) .gt. 0)
				n = n+1
			end do
			newaggregate_hb(n,0) = 2
			newaggregate_hb(n,1) = i
			newaggregate_hb(n,2) = j
		endif
	else 
		call mergeaggregates_hb(indexi, indexj, i, j)
	endif

	return

	end

!*******************************************************************************************

	subroutine mergeaggregates_hb(i, j, chaini, chainj)

	use global

	implicit none

	integer i, j, n, chaini, chainj

	integer sizei, sizej
	
	if (i .eq. 0) then
		newaggregate_hb(j,0) = newaggregate_hb(j,0) + 1
		newaggregate_hb(j, (newaggregate_hb(j,0))) = chaini
		return
	endif

	if (j .eq. 0) then
		newaggregate_hb(i,0) = newaggregate_hb(i,0) + 1
		newaggregate_hb(i, (newaggregate_hb(i,0))) = chainj
		return
	endif

	sizei = newaggregate_hb(i,0)
	sizej = newaggregate_hb(j,0)
	newaggregate_hb(i,0) = sizei + sizej
	newaggregate_hb(j,0) = 0	

	do n = 1, sizej
		newaggregate_hb(i,sizei+n) = newaggregate_hb(j,n)
		newaggregate_hb(j,n) = 0
	end do
		
	return

	end

!**************************************************************************************
