        subroutine write_xyz()

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

	oxygen = chnln1*nc+chnln2*nc2
	write(traj,*) noptotal		!+oxygen
	!txyz = t + tfalse
	write(traj,*) " "
	   
	!do k = 1,(noptotal)
	!	xtmp(k)=sv(1,k)+sv(4,k)*tfalse
	!!	ytmp(k)=sv(2,k)+sv(5,k)*tfalse
	!	ztmp(k)=sv(3,k)+sv(6,k)*tfalse
	!enddo
        	
	do k=1,(noptotal)
         	xtmp(k)=xtmp(k)*boxl
         	ytmp(k)=ytmp(k)*boxl
         	ytmp(k)=ztmp(k)*boxl
      	enddo

	do k=1,nop1/numbeads1+nop2/numbeads2
            if (k .le. nop1/numbeads1) then 
		m =(k-1)*(numbeads1+chnln1)+1
		p = 0
		j =1
		n = (k-1)*numbeads1+1			  
 			do while (j .le. chnln1)								
				write(traj,*) 'N', xtmp(n+chnln1),ytmp(n+chnln1),ztmp(n+chnln1)								
				write(traj,*) 'CA', xtmp(n),ytmp(n),ztmp(n)					
				write(traj,*) 'C', xtmp(n+2*chnln1),ytmp(n+2*chnln1),ztmp(n+2*chnln1)				
   	  			!if(j.ne.chnln1) then
   	  			!	rx1=xtmp(n+2*chnln1)-(xtmp(n)+xtmp(n+1+chnln1))/2.d0
   	  			!	ry1=ytmp(n+2*chnln1)-(ytmp(n)+ytmp(n+1+chnln1))/2.d0
   	  			!	rz1=ztmp(n+2*chnln1)-(ztmp(n)+ztmp(n+1+chnln1))/2.d0
   	  			!else
   	  			!	rx1=xtmp(n+2*chnln1)-xtmp(n)
   	  		        !	ry1=ytmp(n+2*chnln1)-ytmp(n)
   	  			!	rz1=ztmp(n+2*chnln1)-ztmp(n)
   	  			!endif
   	  			!	r1=sqrt(rx1**2+ry1**2+rz1**2)
   	  			!	rx1=xtmp(n+2*chnln1)+1.231d0/r1*rx1
   	  			!	ry1=ytmp(n+2*chnln1)+1.231d0/r1*ry1
   	  			!	rz1=ztmp(n+2*chnln1)+1.231d0/r1*rz1

   	  			!write(traj,*) 'O', rx1, ry1, rz1								
				if (fside1(j).ne.0) then
				  	write(traj,*) 'CB', xtmp(n-p+3*chnln1),ytmp(n-p+3*chnln1),ztmp(n-p+3*chnln1)								  
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
                  		write(traj,*) 'N', xtmp(n+chnln2),ytmp(n+chnln2),ztmp(n+chnln2)								
				write(traj,*) 'CA', xtmp(n),ytmp(n),ztmp(n)					
				write(traj,*) 'C', xtmp(n+2*chnln2),ytmp(n+2*chnln2),ztmp(n+2*chnln2)				
   	  			!if(j.ne.chnln2) then
   	  			!	rx1=xtmp(n+2*chnln2)-(xtmp(n)+xtmp(n+1+chnln2))/2.d0
   	  			!	ry1=ytmp(n+2*chnln2)-(ytmp(n)+ytmp(n+1+chnln2))/2.d0
   	  			!	rz1=ztmp(n+2*chnln2)-(ztmp(n)+ztmp(n+1+chnln2))/2.d0
   	  			!else
   	  			!	rx1=xtmp(n+2*chnln2)-xtmp(n)
   	  		        !	ry1=ytmp(n+2*chnln2)-ytmp(n)
   	  			!	rz1=ztmp(n+2*chnln2)-ztmp(n)
   	  			!endif
   	  			!r1=sqrt(rx1**2+ry1**2+rz1**2)
   	  			!rx1=xtmp(n+2*chnln2)+1.231d0/r1*rx1
   	  			!ry1=ytmp(n+2*chnln2)+1.231d0/r1*ry1
   	  			!rz1=ztmp(n+2*chnln2)+1.231d0/r1*rz1
   	  			!write(traj,*) 'O',rx1,ry1,rz1
				if (fside2(j).ne.0) then
					write(traj,*) 'CB', xtmp(n-p+3*chnln2),ytmp(n-p+3*chnln2),ztmp(n-p+3*chnln2)							  
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