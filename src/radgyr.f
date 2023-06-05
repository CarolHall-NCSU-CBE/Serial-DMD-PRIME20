   	subroutine radgyr(rg_avg)

#include "def.h"

  	use global
	use inputreadin

      	implicit none

      	real*8 comx,comy,comz,rx,ry,rz,rgsum,rxyz(3)
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	real*8 sv_temp(3,(noptotal)),rg_avg
      	integer i,j,k,ii
      	real*8 bkb_rg_sqd((nop1/numbeads1)+(nop2/numbeads2))

!     	calculate each chain's center of mass, then find rad of gyr
      	rg_avg = 0.d0

      	do i=1,numchains
      		if (i .le. nop1/numbeads1) then
	 		ii = (i-1)*numbeads1
         		do j=ii+1,ii+numbeads1
            			sv_temp(1,j)=sv(1,j)+sv(4,j)*tfalse
            			sv_temp(2,j)=sv(2,j)+sv(5,j)*tfalse
            			sv_temp(3,j)=sv(3,j)+sv(6,j)*tfalse
         		end do
	 		do j=ii+2,ii+chnln1
	    			do k = 1, 3
	       			rxyz(k)=sv_temp(k,j)-sv_temp(k,j-1)
	       			if (rxyz(k) .gt. half) then
                  				sv_temp(k,j) = sv_temp(k,j) - boxl
               			else if (rxyz(k) .lt. -1.d0*half) then
                  				sv_temp(k,j) = sv_temp(k,j) + boxl
               			end if
	    			end do
	 		end do
			do j=ii+1,ii+chnln1
	    			do k = 1, 3
	       			rxyz(k)=sv_temp(k,chnln1+j)-sv_temp(k,j)
	       			if (rxyz(k) .gt. half) then
                  				sv_temp(k,chnln1+j) = sv_temp(k,chnln1+j) - boxl
 	       			else if (rxyz(k) .lt. -1.d0*half) then
	          				sv_temp(k,chnln1+j) = sv_temp(k,chnln1+j) + boxl
	       			end if
	    			end do
	 		end do
	 		do j=ii+1,ii+chnln1
	    			do k = 1, 3
               			rxyz(k)=sv_temp(k,2*chnln1+j)-sv_temp(k,j)
	       			if (rxyz(k) .gt. half) then
                  				sv_temp(k,2*chnln1+j) = sv_temp(k,2*chnln1+j) - boxl
               			else if (rxyz(k) .lt. -1.d0*half) then
                  				sv_temp(k,2*chnln1+j) = sv_temp(k,2*chnln1+j) + boxl
               			end if
            			end do
	 		enddo
	         	comx=0.0
         		comy=0.0
         		comz=0.0
         		rgsum=0.0
         		do j=(i-1)*numbeads1+1,(i-1)*numbeads1+3*chnln1
            			comx=comx+sv_temp(1,j)
            			comy=comy+sv_temp(2,j) 
            			comz=comz+sv_temp(3,j) 
         		enddo
	 	       comx=comx/dble(3*chnln1)
         		comy=comy/dble(3*chnln1)
         		comz=comz/dble(3*chnln1)
         		do j=(i-1)*numbeads1+1,(i-1)*numbeads1+3*chnln1
            			rx=(sv_temp(1,j)-comx)**2
            			ry=(sv_temp(2,j)-comy)**2
            			rz=(sv_temp(3,j)-comz)**2
            			rgsum=rgsum+rx+ry+rz
         		enddo
!         		rg_sqd(i)=rgsum/dble(3*chnln1)
	 		rg_avg=rg_avg+rgsum
      		else
	 		ii = nop1+(i-nop1/numbeads1-1)*numbeads2
         		do j=ii+1,ii+numbeads2
            			sv_temp(1,j)=sv(1,j)+sv(4,j)*tfalse
            			sv_temp(2,j)=sv(2,j)+sv(5,j)*tfalse
            			sv_temp(3,j)=sv(3,j)+sv(6,j)*tfalse
         		end do
	 		do j=ii+2,ii+chnln2
	    			do k = 1, 3
	       			rxyz(k)=sv_temp(k,j)-sv_temp(k,j-1)
	       			if (rxyz(k) .gt. half) then
                  				sv_temp(k,j) = sv_temp(k,j) - boxl
               			else if (rxyz(k) .lt. -1.d0*half) then
                  				sv_temp(k,j) = sv_temp(k,j) + boxl
               			end if
	    			end do
	 		end do
         		do j=ii+1,ii+chnln2
	    			do k = 1, 3
	       			rxyz(k)=sv_temp(k,chnln2+j)-sv_temp(k,j)
	       			if (rxyz(k) .gt. half) then
                  				sv_temp(k,chnln2+j) = sv_temp(k,chnln2+j) - boxl
 	       			else if (rxyz(k) .lt. -1.d0*half) then
	          				sv_temp(k,chnln2+j) = sv_temp(k,chnln2+j) + boxl
	       			end if
	    			end do
	 		end do
	 		do j=ii+1,ii+chnln2
	    			do k = 1, 3
               			rxyz(k)=sv_temp(k,2*chnln2+j)-sv_temp(k,j)
	       			if (rxyz(k) .gt. half) then
                  				sv_temp(k,2*chnln2+j) = sv_temp(k,2*chnln2+j) - boxl
               			else if (rxyz(k) .lt. -1.d0*half) then
                  				sv_temp(k,2*chnln2+j) = sv_temp(k,2*chnln2+j) + boxl
               			end if
            			end do
	 		enddo
			comx=0.0
         		comy=0.0
         		comz=0.0
         		rgsum=0.0
         		do j=(i-1)*numbeads2+1,(i-1)*numbeads2+3*chnln2
            			comx=comx+sv_temp(1,j)
            			comy=comy+sv_temp(2,j) 
            			comz=comz+sv_temp(3,j) 
         		enddo
	 		comx=comx/dble(3*chnln2)
         		comy=comy/dble(3*chnln2)
         		comz=comz/dble(3*chnln2)
         		do j=(i-1)*numbeads2+1,(i-1)*numbeads2+3*chnln2
            			rx=(sv_temp(1,j)-comx)**2
            			ry=(sv_temp(2,j)-comy)**2
            			rz=(sv_temp(3,j)-comz)**2
            			rgsum=rgsum+rx+ry+rz
         		enddo
	 		rg_avg=rg_avg+rgsum
      		endif
	enddo

      	rg_avg=dsqrt(rg_avg/dble((nop1/numbeads1)*3*chnln1+(nop2/numbeads2)*3*chnln2))*boxl_orig
      
      	return

      	end
