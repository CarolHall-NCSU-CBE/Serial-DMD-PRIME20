        subroutine write_rasmol

#include "def.h"

        USE GLOBAL
	use inputreadin
                      
	IMPLICIT NONE

        integer i,j,chainn,n,ii,EVCODE
	real*8 RXIJ,RYIJ,RZIJ,RIJSQ,BLMAXSQ,rx1,ry1,rz1,r1
        character*1 col(numchains)
        character*3 amino_name
        integer k2,k4,m
     
7       format(A4,3X,I4,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3)
          k2=0
          k4=0
       !print *, chnln
                col(1)='A'   
                col(2)='B'   
                col(3)='C'   
                col(4)='D'   
                col(5)='E'   
                col(6)='F'   
                col(7)='G'   
                col(8)='H'   
                col(9)='I'   
                col(10)='J'   
                col(11)='K'   
                col(12)='L'   
                col(13)='M'   
                col(14)='N'   
                col(15)='O'   
                col(16)='P'   
                col(17)='Q'   
                col(18)='R'   
                col(19)='S'   
                col(20)='T'   
                col(21)='U'   
                col(22)='V'   
                col(23)='W'   
                col(24)='X'   
   	  	do i=1,nop1,numbeads1
   	  		k2=k2+1
   	  		do j=1,chnln1
				k4=k4+1  
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
   	  	     		write(runpdb,7) 'ATOM',k4,' N  ',amino_name,col(chnnum(i)),j,sv(1,i-1+j+chnln1),sv(2,i-1+j+chnln1),sv(3,i-1+j+chnln1)
                                !print *, 'N',i-1+j+chnln1,sv(1,i-1+j+chnln1),sv(2,i-1+j+chnln1),sv(3,i-1+j+chnln1)
				k4=k4+1
		      		write(runpdb,7) 'ATOM',k4,' CA ',amino_name,col(chnnum(i)),j,sv(1,i-1+j),sv(2,i-1+j),sv(3,i-1+j)
                                !print *, 'C',i-1+j,sv(1,i-1+j),sv(2,i-1+j),sv(3,i-1+j)
				k4=k4+1
		      		write(runpdb,7) 'ATOM',k4,' C  ',amino_name,col(chnnum(i)),j,sv(1,i-1+j+2*chnln1),sv(2,i-1+j+2*chnln1),sv(3,i-1+j+2*chnln1)
   	  			if(j.ne.chnln1) then
   	  				rx1=SV(1,i-1+j+2*chnln1)-(SV(1,i-1+j)+SV(1,i+j+chnln1))/2.d0
   	  				ry1=SV(2,i-1+j+2*chnln1)-(SV(2,i-1+j)+SV(2,i+j+chnln1))/2.d0
   	  				rz1=SV(3,i-1+j+2*chnln1)-(SV(3,i-1+j)+SV(3,i+j+chnln1))/2.d0
   	  			else
   	  				rx1=SV(1,i-1+j+2*chnln1)-SV(1,i-1+j)
   	  		        	ry1=SV(2,i-1+j+2*chnln1)-SV(2,i-1+j)
   	  				rz1=SV(3,i-1+j+2*chnln1)-SV(3,i-1+j)
   	  			endif
   	  			r1=dsqrt(rx1**2+ry1**2+rz1**2)
   	  			rx1=SV(1,i-1+j+2*chnln1)+1.231d0/r1*rx1
   	  			ry1=SV(2,i-1+j+2*chnln1)+1.231d0/r1*ry1
   	  			rz1=SV(3,i-1+j+2*chnln1)+1.231d0/r1*rz1
				k4=k4+1
   	  			write(runpdb,7) 'ATOM',k4,' O  ',amino_name,col(chnnum(i)),j,rx1,ry1,rz1
				if(fside1(j).ne.0) then	
					k4=k4+1
                                       m = floor((j+2)/5.0d0)+floor(j/5.0d0)
		      			write(runpdb,7) 'ATOM',k4,' CB ',amino_name,col(chnnum(i)),j,sv(1,i-1+j+3*chnln1-m),sv(2,i-1+j+3*chnln1-m),sv(3,i-1+j+3*chnln1-m)
				endif
	     		enddo
		enddo
		
   	  	do i=nop1+1,nop1+nop2,numbeads2
   	  		k2=k2+1
   	  		do j=1,chnln2
				k4=k4+1  
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
   	  	     		write(runpdb,7) 'ATOM',k4,' N  ',amino_name,col(chnnum(i)-nop1/chnln1),j,sv(1,i-1+j+chnln2),sv(2,i-1+j+chnln2),sv(3,i-1+j+chnln2)
                                !print *, 'N',i-1+j+chnln2,sv(1,i-1+j+chnln2),sv(2,i-1+j+chnln2),sv(3,i-1+j+chnln2)
				k4=k4+1
		      		write(runpdb,7) 'ATOM',k4,' CA ',amino_name,col(chnnum(i)-nop1/chnln1),j,sv(1,i-1+j),sv(2,i-1+j),sv(3,i-1+j)
                                !print *, 'C',i-1+j,sv(1,i-1+j),sv(2,i-1+j),sv(3,i-1+j)
				k4=k4+1
		      		write(runpdb,7) 'ATOM',k4,' C  ',amino_name,col(chnnum(i)-nop1/chnln1),j,sv(1,i-1+j+2*chnln2),sv(2,i-1+j+2*chnln2),sv(3,i-1+j+2*chnln2)
   	  			if(j.ne.chnln2) then
   	  				rx1=SV(1,i-1+j+2*chnln2)-(SV(1,i-1+j)+SV(1,i+j+chnln2))/2.d0
   	  				ry1=SV(2,i-1+j+2*chnln2)-(SV(2,i-1+j)+SV(2,i+j+chnln2))/2.d0
   	  				rz1=SV(3,i-1+j+2*chnln2)-(SV(3,i-1+j)+SV(3,i+j+chnln2))/2.d0
   	  			else
   	  				rx1=SV(1,i-1+j+2*chnln2)-SV(1,i-1+j)
   	  		        	ry1=SV(2,i-1+j+2*chnln2)-SV(2,i-1+j)
   	  				rz1=SV(3,i-1+j+2*chnln2)-SV(3,i-1+j)
   	  			endif
   	  			r1=dsqrt(rx1**2+ry1**2+rz1**2)
   	  			rx1=SV(1,i-1+j+2*chnln2)+1.231d0/r1*rx1
   	  			ry1=SV(2,i-1+j+2*chnln2)+1.231d0/r1*ry1
   	  			rz1=SV(3,i-1+j+2*chnln2)+1.231d0/r1*rz1
				k4=k4+1
   	  			write(runpdb,7) 'ATOM',k4,' O  ',amino_name,col(chnnum(i)-nop1/chnln1),j,rx1,ry1,rz1
				if(fside2(j).ne.0) then	
					k4=k4+1
                                       m = floor((j+2)/5.0d0)+floor(j/5.0d0)
		      			write(runpdb,7) 'ATOM',k4,' CB ',amino_name,col(chnnum(i)-nop1/chnln1),j,sv(1,i-1+j+3*chnln2-m),sv(2,i-1+j+3*chnln2-m),sv(3,i-1+j+3*chnln2-m)
				endif
	     		enddo
		enddo
		
	!BLMAXSQ = (BDLN(chnln1)*(1+DEL1))**2         !Yiming
	BLMAXSQ = 100
!	do i = 1, NOP-1
!	   do j = i + 1, NOP
!	      EVCODE=EV_CODE(J,I)
!	      IF (EVCODE == 4 .OR. EVCODE == 5 .OR. EVCODE == 6 .OR. EVCODE == 10) THEN
!	         RXIJ=SV(1,I)-SV(1,J)
!                 RYIJ=SV(2,I)-SV(2,J)
!                 RZIJ=SV(3,I)-SV(3,J)
!	         RIJSQ=RXIJ*RXIJ+RYIJ*RYIJ+RZIJ*RZIJ
!	         IF (RIJSQ .LE. BLMAXSQ) THEN
!		    write(runpdb,'(A6,I5,I5)') 'CONECT', i, j
!	         END IF
!	      END IF
!	   end do
!	end do
!           
!	return
           
        end   

