	subroutine write_rasmol

#include "def.h"

     	use global
	use inputreadin
                      
	implicit none

     	integer i,j,chainn,n,ii,evcode
	real*8 rxij,ryij,rzij,rijsq,blmaxsq
     	character*2 col(2)
        
    	do chainn=1, (numchains)-(nop2/numbeads2)
	          	col(1) = 'Al'
              	col(2) = 'K'
           	do i = (chainn-1)*numbeads1+1, chainn*numbeads1
	      		ii = i-((chnnum(i)-1)*numbeads1)
	      		if (identity(i) .gt. 8) then
	         		write(runpdb,'(a6,i5,a3,1x,15x,3f8.3)') 'ATOM  ',i,col(1),sv(1,i),sv(2,i),sv(3,i)
			endif
	      		if (ii==3*chnln1) then
	         		write(runpdb,'(a6,i5,a3,1x,15x,3f8.3)') 'ATOM  ',i,'S',sv(1,i),sv(2,i),sv(3,i)
			endif
			if (identity(i) .le. 8) then
	         		write(runpdb,'(a6,i5,a3,1x,15x,3f8.3)') 'ATOM  ',i,col(2),sv(1,i),sv(2,i),sv(3,i)
	      		end if
	     	end do	
  	end do

    	do chainn=1, (nop2/numbeads2)
	          	col(1) = 'O'
              	col(2) = 'N'
           	do i = nop1+(chainn-1)*numbeads2+1, nop1+chainn*numbeads2
           		ii=i-nop1-((chnnum(i)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
	      		if (identity(i) .gt. 8) then
	         		write(runpdb,'(a6,i5,a3,1x,15x,3f8.3)') 'ATOM  ',i,col(1),sv(1,i),sv(2,i),sv(3,i)
			endif
	      		if (ii==numbeads1+3*chnln2) then
	         		write(runpdb,'(a6,i5,a3,1x,15x,3f8.3)') 'ATOM  ',i,'S',sv(1,i),sv(2,i),sv(3,i)
	      		endif
			if (identity(i) .le. 8) then
	         		write(runpdb,'(a6,i5,a3,1x,15x,3f8.3)') 'ATOM  ',i,col(2),sv(1,i),sv(2,i),sv(3,i)
	      		end if
	     	end do	
    	end do
    	
!LR: Prints out 3rd species to the rasmol file

!	blmaxsq = (bdln(chnln1+1)*(1+del))**2 
	blmaxsq=100

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do i = 1, (noptotal)-1
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
		do j = i + 1, (noptotal)
	 		evcode=ev_code(j,i)
	   		if (evcode == 4 .or. evcode == 5 .or. evcode == 6 .or. evcode == 10) then
	    			rxij=sv(1,i)-sv(1,j)
                 		ryij=sv(2,i)-sv(2,j)
                 		rzij=sv(3,i)-sv(3,j)
	         		rijsq=rxij*rxij+ryij*ryij+rzij*rzij
	         		if (rijsq .le. blmaxsq) then
		    			write(runpdb,'(a6,i5,i5)') 'CONECT', i, j
	         		end if
	      		end if
!			if (j .eq. bptnr(i)) then
!					write(runpdb,'(a6,i5,i5)') 'CONECT', i, j
!			endif
	   	end do
	end do
           
	return
           
   	end   

