        subroutine write_rasmol

#include "def.h"

        USE GLOBAL
	use inputreadin
                      
	IMPLICIT NONE

        !integer i,j,chainn,n,ii,EVCODE
	real*8 RXIJ,RYIJ,RZIJ,RIJSQ,BLMAXSQ,rx1,ry1,rz1,r1
        character*1 col(nop1/numbeads1+nop2/numbeads2)
        character*3 amino_name
        integer i,j,k,m,n,p
		!real*8 newrx(noptotal), newry(noptotal), newrz(noptotal)
     
7       format(A4,3X,I4,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3)
          k = 0
                col(1)='A'   
                col(2)='B'   
           	
	
            do k=1,nop1/numbeads1+nop2/numbeads2
             if (k .le. nop1/numbeads1) then 
			   m =(k-1)*(numbeads1+chnln1)+1
			   p = 0
				j =1
				n = (k-1)*numbeads1+1			  
 					do while (j .le. chnln1)
						
								if(fside1(j).ne.0) then
                                        if(identity(fside1(j)).eq.10) amino_name='ARG'
                                        if(identity(fside1(j)).eq.11) amino_name='ASN'
                                        if(identity(fside1(j)).eq.12) amino_name='ASP'
                                        if(identity(fside1(j)).eq.13) amino_name='GLN'
                                        if(identity(fside1(j)).eq.14) amino_name='GLU'
                                        if(identity(fside1(j)).eq.15) amino_name='HIS'
                                        if(identity(fside1(j)).eq.16) amino_name='LYS'
                                        if(identity(fside1(j)).eq.17) amino_name='PRO'
                                        if(identity(fside1(j)).eq.18) amino_name='SER'
                                        if(identity(fside1(j)).eq.19) amino_name='THR'
                                        if(identity(fside1(j)).eq.20) amino_name='ALA'
                                        if(identity(fside1(j)).eq.21) amino_name='CYS'
                                        if(identity(fside1(j)).eq.22) amino_name='ILE'
                                        if(identity(fside1(j)).eq.23) amino_name='LEU'
                                        if(identity(fside1(j)).eq.24) amino_name='MET'
                                        if(identity(fside1(j)).eq.25) amino_name='PHE'
                                        if(identity(fside1(j)).eq.26) amino_name='TRP'
                                        if(identity(fside1(j)).eq.27) amino_name='TYR'
                                        if(identity(fside1(j)).eq.28) amino_name='VAL'
                                else
                                        amino_name='GLY'
                                endif
								
							  !new2_rx(m)=newrx2(n+chnln2)
								!new2_rx(m+1)=newrx2(n)
								!new2_rx(m+2)=newrx2(n+2*chnln2)
								
                  write(runpdb,7) 'ATOM',m  ,' N  ',amino_name,col(1),j,sv(1,n+chnln1),sv(2,n+chnln1),sv(3,n+chnln1)								
				  write(runpdb,7) 'ATOM',m+1,' CA ',amino_name,col(1),j,sv(1,n),sv(2,n),sv(3,n)					
				  write(runpdb,7) 'ATOM',m+2,' C  ',amino_name,col(1),j,sv(1,n+2*chnln1),sv(2,n+2*chnln1),sv(3,n+2*chnln1)				
   	  			if(j.ne.chnln1) then
   	  				rx1=SV(1,n+2*chnln1)-(SV(1,n)+SV(1,n+1+chnln1))/2.d0
   	  				ry1=SV(2,n+2*chnln1)-(SV(2,n)+SV(2,n+1+chnln1))/2.d0
   	  				rz1=SV(3,n+2*chnln1)-(SV(3,n)+SV(3,n+1+chnln1))/2.d0
   	  			else
   	  				rx1=SV(1,n+2*chnln1)-SV(1,n)
   	  		        ry1=SV(2,n+2*chnln1)-SV(2,n)
   	  				rz1=SV(3,n+2*chnln1)-SV(3,n)
   	  			endif
   	  			r1=dsqrt(rx1**2+ry1**2+rz1**2)
   	  			rx1=SV(1,n+2*chnln1)+1.231d0/r1*rx1
   	  			ry1=SV(2,n+2*chnln1)+1.231d0/r1*ry1
   	  			rz1=SV(3,n+2*chnln1)+1.231d0/r1*rz1

   	  			  write(runpdb,7) 'ATOM',m+3,' O  ',amino_name,col(1),j,rx1,ry1,rz1								
							if (fside1(j).ne.0) then
								!new2_rx(m+3)=newrx2(n-p+3*chnln2)
				  write(runpdb,7) 'ATOM',m+4,' CB ',amino_name,col(1),j,sv(1,n-p+3*chnln1),sv(2,n-p+3*chnln1),sv(3,n-p+3*chnln1)
															  
							    m = m + 5
						    else 
								p = p + 1							
								m = m + 4
						    endif							
							j = j + 1
							n = n + 1
						end do  					                        

			  
		      else
	
		    !do k=1,nop2/numbeads2
			   m =nop1+chnln1*(nop1/numbeads1)+(k-nop1/numbeads1-1)*(numbeads2+chnln2)+1
			   p = 0
				j =1
				n = nop1+(k-nop1/numbeads1-1)*(numbeads2)+1
						do while (j .le. chnln2)
						
								if(fside2(j).ne.0) then
                                        if(identity(fside2(j)).eq.10) amino_name='ARG'
                                        if(identity(fside2(j)).eq.11) amino_name='ASN'
                                        if(identity(fside2(j)).eq.12) amino_name='ASP'
                                        if(identity(fside2(j)).eq.13) amino_name='GLN'
                                        if(identity(fside2(j)).eq.14) amino_name='GLU'
                                        if(identity(fside2(j)).eq.15) amino_name='HIS'
                                        if(identity(fside2(j)).eq.16) amino_name='LYS'
                                        if(identity(fside2(j)).eq.17) amino_name='PRO'
                                        if(identity(fside2(j)).eq.18) amino_name='SER'
                                        if(identity(fside2(j)).eq.19) amino_name='THR'
                                        if(identity(fside2(j)).eq.20) amino_name='ALA'
                                        if(identity(fside2(j)).eq.21) amino_name='CYS'
                                        if(identity(fside2(j)).eq.22) amino_name='ILE'
                                        if(identity(fside2(j)).eq.23) amino_name='LEU'
                                        if(identity(fside2(j)).eq.24) amino_name='MET'
                                        if(identity(fside2(j)).eq.25) amino_name='PHE'
                                        if(identity(fside2(j)).eq.26) amino_name='TRP'
                                        if(identity(fside2(j)).eq.27) amino_name='TYR'
                                        if(identity(fside2(j)).eq.28) amino_name='VAL'
                                else
                                        amino_name='GLY'
                                endif
								
							  !new2_rx(m)=newrx2(n+chnln2)
								!new2_rx(m+1)=newrx2(n)
								!new2_rx(m+2)=newrx2(n+2*chnln2)
								
                  write(runpdb,7) 'ATOM',m  ,' N  ',amino_name,col(2),j,sv(1,n+chnln2),sv(2,n+chnln2),sv(3,n+chnln2)								
				  write(runpdb,7) 'ATOM',m+1,' CA ',amino_name,col(2),j,sv(1,n),sv(2,n),sv(3,n)					
				  write(runpdb,7) 'ATOM',m+2,' C  ',amino_name,col(2),j,sv(1,n+2*chnln2),sv(2,n+2*chnln2),sv(3,n+2*chnln2)				
   	  			if(j.ne.chnln2) then
   	  				rx1=SV(1,n+2*chnln2)-(SV(1,n)+SV(1,n+1+chnln2))/2.d0
   	  				ry1=SV(2,n+2*chnln2)-(SV(2,n)+SV(2,n+1+chnln2))/2.d0
   	  				rz1=SV(3,n+2*chnln2)-(SV(3,n)+SV(3,n+1+chnln2))/2.d0
   	  			else
   	  				rx1=SV(1,n+2*chnln2)-SV(1,n)
   	  		        ry1=SV(2,n+2*chnln2)-SV(2,n)
   	  				rz1=SV(3,n+2*chnln2)-SV(3,n)
   	  			endif
   	  			r1=dsqrt(rx1**2+ry1**2+rz1**2)
   	  			rx1=SV(1,n+2*chnln2)+1.231d0/r1*rx1
   	  			ry1=SV(2,n+2*chnln2)+1.231d0/r1*ry1
   	  			rz1=SV(3,n+2*chnln2)+1.231d0/r1*rz1

   	  			  write(runpdb,7) 'ATOM',m+3,' O  ',amino_name,col(2),j,rx1,ry1,rz1

				  
							if (fside2(j).ne.0) then
								!new2_rx(m+3)=newrx2(n-p+3*chnln2)
				  write(runpdb,7) 'ATOM',m+4,' CB ',amino_name,col(2),j,sv(1,n-p+3*chnln2),sv(2,n-p+3*chnln2),sv(3,n-p+3*chnln2)
															  
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
	end		