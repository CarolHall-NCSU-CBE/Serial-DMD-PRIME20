        subroutine write_traj()

#include "def.h"

        USE GLOBAL
	use inputreadin
                      
	IMPLICIT NONE

        !integer i,j,chainn,n,ii,EVCODE
	real*4 RXIJ,RYIJ,RZIJ,RIJSQ,BLMAXSQ,rx1,ry1,rz1,r1
	real*4 xtmp(noptotal),ytmp(noptotal),ztmp(noptotal)
        character*1 col(nop1/numbeads1+nop2/numbeads2)
        character*3 amino_name
        integer i,j,k,m,n,p

	integer :: oxygen
	real*4 :: txyz

		   
	do k = 1,(noptotal)
		xtmp(k)=sv(1,k)+sv(4,k)*tfalse
		ytmp(k)=sv(2,k)+sv(5,k)*tfalse
		ztmp(k)=sv(3,k)+sv(6,k)*tfalse
	enddo
        	
	do k=1,(noptotal)
         	xtmp(k)=xtmp(k)*boxl_orig
         	ytmp(k)=ytmp(k)*boxl_orig
         	ztmp(k)=ztmp(k)*boxl_orig
      	enddo

	do k=1,nop1/numbeads1+nop2/numbeads2
            if (k .le. nop1/numbeads1) then 
		m =(k-1)*(numbeads1+chnln1)+1
		p = 0
		j =1
		n = (k-1)*numbeads1+1			  
 			do while (j .le. chnln1)								
				write(traj) xtmp(n+chnln1),ytmp(n+chnln1),ztmp(n+chnln1)								
				write(traj) xtmp(n),ytmp(n),ztmp(n)					
				write(traj) xtmp(n+2*chnln1),ytmp(n+2*chnln1),ztmp(n+2*chnln1)				
   	  											
				if (fside1(j).ne.0) then
				  	write(traj) xtmp(n-p+3*chnln1),ytmp(n-p+3*chnln1),ztmp(n-p+3*chnln1)								  
					m = m + 5
				else 
					p = p + 1							
					m = m + 4
				endif							
					j = j + 1
					n = n + 1
			end do  					                          
		else
			m =nop1+chnln1*(nop1/numbeads1)+(k-nop1/numbeads1-1)*(numbeads2+chnln2)+1
			p = 0
			j =1
			n = nop1+(k-nop1/numbeads1-1)*(numbeads2)+1
			do while (j .le. chnln2)																
                  		write(traj) xtmp(n+chnln2),ytmp(n+chnln2),ztmp(n+chnln2)								
				write(traj) xtmp(n),ytmp(n),ztmp(n)					
				write(traj) xtmp(n+2*chnln2),ytmp(n+2*chnln2),ztmp(n+2*chnln2)				
   				if (fside2(j).ne.0) then
					write(traj) xtmp(n-p+3*chnln2),ytmp(n-p+3*chnln2),ztmp(n-p+3*chnln2)							  
					m = m + 5
				else 
					p = p + 1							
					m = m + 4
				endif							
					j = j + 1
					n = n + 1
			end do  					
                end if       				
	end do
	return
	end		