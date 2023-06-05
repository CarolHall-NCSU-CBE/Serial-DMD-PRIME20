	subroutine calcaggregates_hh(num_aggs_hh)

#include "def.h"

	use global
	use inputreadin

	implicit none

	integer i, j, n, m, p, num_aggs_hh

	do i = 1, numchains/2
		do j = 0, numchains
			newaggregate_hh(i,j) = 0
		end do
	end do

	do i = 1, numchains - 1
		do j = i+1, numchains
			if (matrix_hh(i,j) .ne. 0) then
				call addtoaggregate_hh(i,j) 
			endif
		enddo
	enddo

	num_aggs_hh = 0

	do n = 1, numchains/2
		if (newaggregate_hh(n,0) .gt. 0) then
			num_aggs_hh = num_aggs_hh + 1
			aggregate_hh(num_aggs_hh,0)=newaggregate_hh(n,0)
			do m = 1, newaggregate_hh(n,0)
				aggregate_hh(num_aggs_hh,m) = newaggregate_hh(n,m)
			end do
		endif
	end do

	return

	end
		    
!*********************************************************************************************
		
	subroutine addtoaggregate_hh(i,j)

	use global

	implicit none

	integer i, j, n, m, indexi, indexj

	indexi = 0
	indexj = 0

	do n = 1, numchains/2
		do m = 1, newaggregate_hh(n,0)
			if (newaggregate_hh(n,m) .eq. i) then
				indexi = n
			else if (newaggregate_hh(n,m) .eq. j) then
				indexj = n
			endif
		end do
	end do

	if (indexi .eq. indexj) then
		if (indexi .eq. 0) then
			n = 1
			do while (newaggregate_hh(n,0) .gt. 0)
				n = n+1
			end do
			newaggregate_hh(n,0) = 2
			newaggregate_hh(n,1) = i
			newaggregate_hh(n,2) = j
		endif
	else 
		call mergeaggregates_hh(indexi, indexj, i, j)
	endif

	return

	end

!*******************************************************************************************

	subroutine mergeaggregates_hh(i, j, chaini, chainj)

	use global

	implicit none

	integer i, j, n, chaini, chainj
	integer sizei, sizej
	
	if (i .eq. 0) then
		newaggregate_hh(j,0) = newaggregate_hh(j,0) + 1
		newaggregate_hh(j, (newaggregate_hh(j,0))) = chaini
		return
	endif

	if (j .eq. 0) then
		newaggregate_hh(i,0) = newaggregate_hh(i,0) + 1
		newaggregate_hh(i, (newaggregate_hh(i,0))) = chainj
		return
	endif

	sizei = newaggregate_hh(i,0)
	sizej = newaggregate_hh(j,0)
	newaggregate_hh(i,0) = sizei + sizej
	newaggregate_hh(j,0) = 0	

	do n = 1, sizej
		newaggregate_hh(i,sizei+n) = newaggregate_hh(j,n)
		newaggregate_hh(j,n) = 0
	end do
		
	return

	end

!**************************************************************************************
