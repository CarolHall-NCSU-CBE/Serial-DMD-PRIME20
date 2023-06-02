	subroutine cell_add

#include "def.h"

	use global
	use inputreadin

	implicit none

	integer k, x, y, z, cell_k, i,j,m
	
	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do k = 1, (noptotal)	
		clinks(k)=0
	end do

	do k = 1, num_cell**3+1
		cell(k)=0
	end do

	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do k = 1, (noptotal)
		x = int((sv(1,k) + half)/width) + n_wrap
		y = int((sv(2,k) + half)/width) + n_wrap
		z = int((sv(3,k) + half)/width) + n_wrap
		cell_k = 1 + x + y*num_cell + z*num_cell*num_cell
		clinks(k)=cell(cell_k)
		cell(cell_k) = k
	end do         

#ifdef debugging

	m = 0

	do k=1, num_cell**3+1
		i=cell(k)
		j = 0
		do while (i .ne. 0)
			j = j + 1
			i = clinks(i)
		end do
		m = m + j
	end do
	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	if (m .ne. (noptotal)) print*, 'total', m

#endif

	end
