	subroutine order_param(num_aggs_hb,num_aggs_hh,ss,ss_hb,ss_hh,parallel)

#include "def.h"

      	use global
	use inputreadin

      	implicit none

      	integer i,j,ii,jj,k,n,nn,m,b(4),n_vector,nn_vector,num_aggs_hb,num_aggs_hh,n_angle
      	real*8 x(numchains,4),y(numchains,4),z(numchains,4),xx,yy,zz,angle,poly
      	real*8 s(numchains),ss,ss_hb,ss_hh,parallel,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,a2,b2,c2,poly2,sx,sy,sz

      	do k = 1, numchains
	 	b(1) = (k-1)*numbeads + 1
	 	b(2) = (k-1)*numbeads + chnln/3 + 1
	 	b(3) = (k-1)*numbeads + (chnln/3)*2 + 1
	 	b(4) = (k-1)*numbeads + chnln
		do m = 1, 4
	    		i = b(m)
	    		x(k,m)=sv(1,i) + sv(4,i)*tfalse
	    		y(k,m)=sv(2,i) + sv(5,i)*tfalse
	    		z(k,m)=sv(3,i) + sv(6,i)*tfalse
	 	end do
      	end do

    	ss_hb=0.d0

   	if (num_aggs_hb .ne. 0) then
         	n_angle = 0
         	nn_vector = 0
         	do nn = 1, num_aggs_hb	    
            		n_vector = 0
            		s(nn)=0.d0
            		do ii=1,aggregate_hb(nn,0)
	       		i = aggregate_hb(nn,ii)
	       		do m = 1, 3
                  			do jj=ii+1,aggregate_hb(nn,0)
	             				j=aggregate_hb(nn,jj)
	             				do n = 1, 3
	                				x1=x(i,m)
                       				y1=y(i,m)
                        				z1=z(i,m)
                        				x2=x(i,m+1)
                        				y2=y(i,m+1)
                        				z2=z(i,m+1)
                        				sx=x(j,n)-x(i,m)
                        				sy=y(j,n)-y(i,m)
                        				sz=z(j,n)-z(i,m)
	                				sx=sx-dnint(sx)
	                				sy=sy-dnint(sy)
	                				sz=sz-dnint(sz)
                        				x3=x(j,n)-sx
                        				y3=y(j,n)-sy
                        				z3=z(j,n)-sz
	                				x3=x3-dnint(x3)
	                				y3=y3-dnint(y3)
	                				z3=z3-dnint(z3)
                        				x4=x(j,n+1)-sx
                        				y4=y(j,n+1)-sy
                        				z4=z(j,n+1)-sz
	                				x4=x4-dnint(x4)
	                				y4=y4-dnint(y4)
	                				z4=z4-dnint(z4)
		        				xx=x2-x1
		        				yy=y2-y1
		        				zz=z2-z1
		        				xx=xx-dnint(xx)
		        				yy=yy-dnint(yy)
		        				zz=zz-dnint(zz)
                        				a2=xx**2+yy**2+zz**2
		        				xx=x4-x3
		        				yy=y4-y3
		        				zz=z4-z3
		        				xx=xx-dnint(xx)
		        				yy=yy-dnint(yy)
		        				zz=zz-dnint(zz)
                        				b2=xx**2+yy**2+zz**2
	                				xx=x4-x2
                        				yy=y4-y2
                        				zz=z4-z2
                        				xx=xx-dnint(xx)
                        				yy=yy-dnint(yy)
                        				zz=zz-dnint(zz)
                        				c2=xx**2+yy**2+zz**2
                        				poly=(a2+b2-c2)/(2.d0*dsqrt(a2)*dsqrt(b2))
	                				if (poly .gt. 1.d0) poly = 1.d0
	                				if (poly .lt. -1.d0) poly = -1.d0
                        				poly2=3.d0*poly**2
                        				angle=acos(poly)
	                				angle=180.d0*angle/pi
                        				if (angle .lt. 90.d0) then 
	                   					n_angle=n_angle+1
	                				end if
		        				nn_vector=nn_vector+1
                        				s(nn)=s(nn)+(poly2-1.d0)/2.d0
                        				n_vector=n_vector+1
	             				end do	   	
	          			end do
               		end do  
            		end do
	    		if (n_vector .ne. 0) s(nn)=s(nn)/dble(n_vector)
	    		ss_hb = ss_hb + s(nn)
         	end do
         	ss_hb=ss_hb/dble(num_aggs_hb)  
	 	parallel=dble(n_angle)/nn_vector
      	end if

      	ss_hh=0.d0

      	if (num_aggs_hh .ne. 0) then
         	do nn = 1, num_aggs_hh	    
            		n_vector = 0
            		s(nn)=0.d0
            		do ii=1,aggregate_hh(nn,0)
	       		i = aggregate_hh(nn,ii)
	       		do m = 1, 3
                  			do jj=ii+1,aggregate_hh(nn,0)
	             				j=aggregate_hh(nn,jj)
	             				do n = 1, 3
	                				x1=x(i,m)
                       				y1=y(i,m)
                        				z1=z(i,m)
                        				x2=x(i,m+1)
                        				y2=y(i,m+1)
                        				z2=z(i,m+1)
                        				sx=x(j,n)-x(i,m)
                        				sy=y(j,n)-y(i,m)
                        				sz=z(j,n)-z(i,m)
	                				sx=sx-dnint(sx)
	                				sy=sy-dnint(sy)
	                				sz=sz-dnint(sz)
                        				x3=x(j,n)-sx
                        				y3=y(j,n)-sy
                        				z3=z(j,n)-sz
	                				x3=x3-dnint(x3)
	                				y3=y3-dnint(y3)
	                				z3=z3-dnint(z3)
                        				x4=x(j,n+1)-sx
                        				y4=y(j,n+1)-sy
                        				z4=z(j,n+1)-sz
	                				x4=x4-dnint(x4)
	                				y4=y4-dnint(y4)
	                				z4=z4-dnint(z4)
		        				xx=x2-x1
		        				yy=y2-y1
		        				zz=z2-z1
		        				xx=xx-dnint(xx)
		        				yy=yy-dnint(yy)
		        				zz=zz-dnint(zz)
                        				a2=xx**2+yy**2+zz**2
		        				xx=x4-x3
		        				yy=y4-y3
		        				zz=z4-z3
		        				xx=xx-dnint(xx)
		        				yy=yy-dnint(yy)
		        				zz=zz-dnint(zz)
                        				b2=xx**2+yy**2+zz**2
	                				xx=x4-x2
                        				yy=y4-y2
                        				zz=z4-z2
                        				xx=xx-dnint(xx)
                        				yy=yy-dnint(yy)
                        				zz=zz-dnint(zz)
                        				c2=xx**2+yy**2+zz**2
                        				poly=(a2+b2-c2)/(2.d0*dsqrt(a2)*dsqrt(b2))
	                				if (poly .gt. 1.d0) poly = 1.d0
	                				if (poly .lt. -1.d0) poly = -1.d0
                        				poly2=3.d0*poly**2
                        				s(nn)=s(nn)+(poly2-1.d0)/2.d0
                        				n_vector=n_vector+1
	             				end do	   	
	          			end do
               		end do  
            		end do
	    		if (n_vector .ne. 0) s(nn)=s(nn)/dble(n_vector)
	    		ss_hh = ss_hh + s(nn)
         	end do
         	ss_hh=ss_hh/dble(num_aggs_hh)  
      	end if

      	n_vector = 0
      	ss=0.d0

      	do i=1,numchains-1
         	do m = 1, 3
            		do j=i+1,numchains
               		do n = 1, 3
                  			x1=x(i,m)
                  			y1=y(i,m)
                  			z1=z(i,m)
                  			x2=x(i,m+1)
                  			y2=y(i,m+1)
                  			z2=z(i,m+1)
                  			sx=x(j,n)-x(i,m)
                  			sy=y(j,n)-y(i,m)
                  			sz=z(j,n)-z(i,m)
	          			sx=sx-dnint(sx)
	          			sy=sy-dnint(sy)
	          			sz=sz-dnint(sz)
                  			x3=x(j,n)-sx
                  			y3=y(j,n)-sy
                  			z3=z(j,n)-sz
	          			x3=x3-dnint(x3)
	          			y3=y3-dnint(y3)
	          			z3=z3-dnint(z3)
                  			x4=x(j,n+1)-sx
                  			y4=y(j,n+1)-sy
                  			z4=z(j,n+1)-sz
	          			x4=x4-dnint(x4)
	          			y4=y4-dnint(y4)
	          			z4=z4-dnint(z4)
		  			xx=x2-x1
		  			yy=y2-y1
		  			zz=z2-z1
		  			xx=xx-dnint(xx)
		  			yy=yy-dnint(yy)
		  			zz=zz-dnint(zz)
                  			a2=xx**2+yy**2+zz**2
		  			xx=x4-x3
		  			yy=y4-y3
		  			zz=z4-z3
		  			xx=xx-dnint(xx)
		  			yy=yy-dnint(yy)
		  			zz=zz-dnint(zz)
                  			b2=xx**2+yy**2+zz**2
	          			xx=x4-x2
                  			yy=y4-y2
                  			zz=z4-z2
                  			xx=xx-dnint(xx)
                  			yy=yy-dnint(yy)
                  			zz=zz-dnint(zz)
                  			c2=xx**2+yy**2+zz**2
                  			poly2=3.d0*((a2+b2-c2)/(2.d0*dsqrt(a2)*dsqrt(b2)))**2
                  			ss=ss+(poly2-1.d0)/2.d0
                  			n_vector=n_vector+1
               		end do
            		end do
         	end do
      	end do

      	if (n_vector .ne. 0) ss=ss/dble(n_vector)

      	return
      
      	end
