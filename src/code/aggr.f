	subroutine calcaggregates(coll)

#include "def.h"

	use global
	use inputreadin

	implicit none

	integer*8 coll
	integer i, j, n, m, p, num_aggs,agg_chains(numchains),aggregated_nchains
	real agg_size(numchains)

	do i = 1, numchains
		agg_size(i) = 0.0                
	end do

	do i = 1, numchains/2
		do j = 0, numchains
			newaggregate(i,j) = 0
		end do
	end do

	do i = 1, numchains - 1
		do j = i+1, numchains
			if ((matrix_hb(i,j) .ne. 0) .or. (matrix_hh(i,j) .ne. 0)) then
				call addtoaggregate(i,j) 
			endif
		enddo
	enddo

	num_aggs = 0

	do n = 1, numchains/2
		if (newaggregate(n,0) .gt. 0) then
			num_aggs = num_aggs + 1
			aggregate(num_aggs,0)=newaggregate(n,0)
			do m = 1, newaggregate(n,0)
				aggregate(num_aggs,m) = newaggregate(n,m)
			end do
		endif
	end do

	do n = 1, numchains/2
		if (newaggregate(n,0) .gt. 0) then
			agg_size(newaggregate(n,0)) = agg_size(newaggregate(n,0)) + newaggregate(n,0) 
		end if 
	end do

	do n = 1, numchains
		if (num_aggs .ne. 0) agg_size(n) = agg_size(n)/dble(numchains)*100.0	   
	end do

	do i = 1, numchains
		agg_chains(i) = 0
	end do  
      
	do i = 1, numchains
		do j = i+1, numchains
			if ((matrix_hb(i,j) .ne. 0) .or. (matrix_hh(i,j) .ne. 0)) then
				agg_chains(i) = 1
				agg_chains(j) = 1
			endif   
		end do
	end do

	aggregated_nchains = 0

	do i = 1, numchains
		aggregated_nchains = aggregated_nchains + agg_chains(i)
	end do

	write(803,100) coll,t,dble(aggregated_nchains)/dble(aggregated_nchains)*100.0,agg_size
100    format(i15,2f12.4, 1000f7.2)

	return

	end
		    
!*********************************************************************************************
		
	subroutine addtoaggregate(i,j)

	use global

	integer i, j, n, m, indexi, indexj

	indexi = 0
	indexj = 0

	do n = 1, numchains/2
		do m = 1, newaggregate(n,0)
			if (newaggregate(n,m) .eq. i) then
				indexi = n
			else if (newaggregate(n,m) .eq. j) then
				indexj = n
			endif
		end do
	end do

	if (indexi .eq. indexj) then
		if (indexi .eq. 0) then
			n = 1
			do while (newaggregate(n,0) .gt. 0)
				n = n+1
			end do
			newaggregate(n,0) = 2
			newaggregate(n,1) = i
			newaggregate(n,2) = j
		endif
	else 
		call mergeaggregates(indexi, indexj, i, j)
	endif

	return

	end

!*******************************************************************************************

	subroutine mergeaggregates(i, j, chaini, chainj)

	use global

	integer i, j, n, chaini, chainj
	integer sizei, sizej
	
	if (i .eq. 0) then
		newaggregate(j,0) = newaggregate(j,0) + 1
		newaggregate(j, (newaggregate(j,0))) = chaini
		return
	endif

	if (j .eq. 0) then
		newaggregate(i,0) = newaggregate(i,0) + 1
		newaggregate(i, (newaggregate(i,0))) = chainj
		return
	endif

	sizei = newaggregate (i,0)
	sizej = newaggregate (j,0)
	newaggregate(i,0) = sizei + sizej
	newaggregate(j,0) = 0	

	do n = 1, sizej
		newaggregate(i,sizei+n) = newaggregate(j,n)
		newaggregate(j,n) = 0
	end do
		
	return

	end

!**************************************************************************************


