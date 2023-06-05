!	subroutine nbor determines and stores neighbors

	subroutine nbor()

#include "def.h"

       use global
	use inputreadin

	implicit none

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	integer evcode,i,j,k,l,c,bead,i_n,j_n,n_bead,f_bead,n_cell(noptotal),n,ncnbor
	real*8 rxij,ryij,rzij,rijsq

!      npt_dn is initialized at start of main and does not change.  
!      npt_dn(x) indicates the position of the first down-neighbor of 
!      x in the array called dnnab (unused spaces are zeroes)

!      nnabdn(x) points to the next available space in dnnab in which 
!      a down-neighbor may be added;  it'll change during dnnab growth

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do l=1,(noptotal)
		na_npt(l)=npt(l)
!LR: Cleaned up neighbor-assignment code
		nnabdn(l)=npt_dn(l)
	enddo
	
!write(fileout,*)'update!',coll
!	calculate distance between each pair of mcs and check to see if
!	each j mc is within the cutoff distance for i mc 
	
	call cell_add

	do c=n_wrap*num_cell*num_cell, num_cell*num_cell*(num_cell-n_wrap)
		n_bead=0
     		bead=cell(c)
           	if (bead .eq. 0) cycle
	   	do while (bead .ne. 0)
	      		n_bead=n_bead+1
	      		n_cell(n_bead)=bead
	      		bead=clinks(bead)
	   	end do
	   	f_bead=n_bead
           	do n = 1, n_nab_cell
	      		ncnbor = c + map(n)
	      		ncnbor = wrap_map(ncnbor)	      
	      		bead=cell(ncnbor)
              	do while (bead .ne. 0)
                 		n_bead=n_bead+1 
                 		n_cell(n_bead)=bead
                 		bead=clinks(bead)
	      		enddo
           	end do
		do i_n=1,f_bead
	      		i = n_cell(i_n)
	      		do j_n=i_n+1,n_bead
	         		j = n_cell(j_n)
				evcode=ev_code(j,i)	
				if (((evcode.ge.4).and.(evcode.le.12)).or.((evcode.ge.17).and.(evcode.lt.27))) then
#ifdef glycine 
                    			if ((bdln_dummy(i) .lt. 99.0).and.(bdln_dummy(j).lt.99.0)) then
#endif              
				
				

                       		if ((j .gt. i)) then
                          			l=na_npt(i)
						if(l .lt. i*maxnbs) then
                          			nb(l)=j
                          			na_npt(i)=l+1
						endif
                          			l=nnabdn(j) 
						if(l .lt. j*maxnbs) then 
                          			dnnab(l)=i   
                          			nnabdn(j)=l+1
						endif
						
                       		else
		          			l=na_npt(j)
						if(l .lt. j*maxnbs) then 
                          			nb(l)=i
                          			na_npt(j)=l+1
						endif
                          			l=nnabdn(i)
						if(l .lt. i*maxnbs) then
                          			dnnab(l)=j
                          			nnabdn(i)=l+1
						endif
                       		endif
				

#ifdef glycine        
                    			endif
#endif                 
                 		else 
	       	    		rxij=sv(1,i)-sv(1,j)
                    			ryij=sv(2,i)-sv(2,j)
                    			rzij=sv(3,i)-sv(3,j)
                    			rxij=rxij-dnint(rxij)
                    			ryij=ryij-dnint(ryij)
                    			rzij=rzij-dnint(rzij)
                    			rijsq=rxij*rxij+ryij*ryij+rzij*rzij
!                   			if j is within cutoff, save it as a neighbor (in array nb)
                    			if (rijsq .le. rlsq(evcode)) then
#ifdef glycine 	      
                       			if ((bdln_dummy(i) .lt. 99.0).and.(bdln_dummy(j).lt.99.0)) then
#endif
		          			if (j .gt. i) then
                             			l=na_npt(i)
                             			nb(l)=j
                             			na_npt(i)=l+1
                             			l=nnabdn(j)
                             			dnnab(l)=i
                             			nnabdn(j)=l+1
                          			else
                             			l=na_npt(j)
                             			nb(l)=i
                             			na_npt(j)=l+1
                             			l=nnabdn(i)
                             			dnnab(l)=j   
                             			nnabdn(i)=l+1
                          			endif
#ifdef glycine
                       			endif
#endif
					endif
				end if
			enddo
		enddo
	enddo

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do k=1,(noptotal)
		na_npt(k)=na_npt(k)-npt(k)
		nnabdn(k)=nnabdn(k)-npt_dn(k)
	enddo

	return

	end



