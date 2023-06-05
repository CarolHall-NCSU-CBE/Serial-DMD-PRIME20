	subroutine cell_link

#include "def.h"

	use global
	use inputreadin

	implicit none

	integer i_cell,x,y,z,ix,iy,iz,d_cell,dx,dy,dz,min_ix,i_map
	integer, dimension(5) :: o = (/ 0, 1, -1, 2, -2 /)
	integer, dimension(3) :: oy = (/ 0, 1, 2 /)

	i_cell(ix,iy,iz) = 1 + mod(ix-1+num_cell,num_cell)+mod(iy-1+num_cell,num_cell)*num_cell &
&                         + mod(iz-1+num_cell,num_cell)*num_cell*num_cell
		
	d_cell = 1

	do iz = 1, 3
		dz=o(iz)
		do iy = 1, 2
			dy=o(iy)
			do ix = 1, 3
				dx=o(ix)
				if (dy==0) then
					if (dz < 1) then
						if (dx < 1) cycle  
                    			else if (dx < 0) then
						cycle
					end if
				end if
				map(d_cell) = dx+(dy+dz*num_cell)*num_cell
				d_cell = d_cell + 1
			end do
		end do
	end do

#if n_wrap==2

	do iz = 1, 5
		dz=o(iz)
		do iy = 1, 3
			dy=oy(iy)
			if (iz<4 .and. iy/=3) then
				min_ix=4
			else
				min_ix=1
			end if   
			do ix = min_ix, 5
				dx=o(ix)
				if (dy==0) then
					if (dz < 1) then
						if (dx < 1) cycle
					else if (dx < 0) then   
						cycle
					end if
				end if
				map(d_cell) = dx+(dy+dz*num_cell)*num_cell
				d_cell = d_cell + 1
			end do
		end do
	end do

#endif

	do ix = 1,num_cell
		if (ix <= n_wrap) then
			x = ix + num_cell - 2*n_wrap
		else if (num_cell-ix < n_wrap) then
			x = ix - num_cell + 2*n_wrap
		else   
			x = ix
		end if
		do iy = 1, num_cell
			if (iy <= n_wrap) then
                 		y = iy + num_cell - 2*n_wrap
			else if (num_cell-iy < n_wrap) then  
                 		y = iy - num_cell + 2*n_wrap
              	else
                 		y = iy  
              	end if     
              	do iz = 1, num_cell
                 		if (iz <= n_wrap) then
                    			z = iz + num_cell - 2*n_wrap
                 		else if (num_cell-iz < n_wrap) then
                   			z = iz - num_cell + 2*n_wrap
                 		else
                    			z = iz
                 		end if
                 		i_map = i_cell(ix,iy,iz) 
                 		d_cell = i_cell(x,y,z)
                 		wrap_map(i_map)=d_cell
			end do
		end do
	end do

	end

