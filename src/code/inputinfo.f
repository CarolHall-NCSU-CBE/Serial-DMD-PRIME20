!	subroutine inputinfo accesses the positions and velocities of the 
!	beads as well as the arrays identity, aa, bdln, hp, chnnum,
!	sigma, welldia, epsilon  

	subroutine inputinfo()
   
#include "def.h"

	use global
	use inputreadin
 
	implicit none
	
!LR: Added a third species numbeads variable. I did not have a numbeadstotal, due to those variables only being used rarely in the code
	integer nn,k,l,aa(numbeads1+numbeads2)
      	integer i,iii,iiii,jjjj !! by mookyung
      	real*8 sumvel,tred,const,sumx,sumy,sumz,x1,x2,y1,y2
      	real*8 eptemp,bdtemp,wltemp,bmass_temp !! by mookyung
      	real*8 drca(20),drnh(20),drco(20) !! by mookyung
      	real*8 del_rca(20),del_rnh(20),del_rco(20) !! by mookyung
      	real*8 scalerca,rx2,ry2,rz2,rca2,rnh2,rco2,delx  !! by mookyung
      	real*8 bmass(28)
      	real*8 sz6,sz7,sz8,sz9,sz10
      	character*64 filename

!     	read in starting positions and velocities
      	write(6,*)' '

#ifndef runr
	write(6,*)'starting from an ideal conformation'
#ifdef runh
!     	accessing ideal helix as starting configuration 
      	open(unit=7,file='parameters/helix.inp',status='unknown')
      	write(6,*)'starting in parameters/helix.inp conformation'
#elif runb
!     	accessing ideal helix as starting configuration (eight 16mer)
!      	open(unit=7,file='parameters/bundles.inp',status='unknown')
      	write(6,*)'starting in parameters/bundles.inp conformation'
#endif

	read (parabundles,*) nn
      	read (parabundles,*) boxl
      	
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	if (nn.ne.(noptotal)) then
		write(6,*)'nop error in reading starting pos and vel '
		!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
         	write(6,*)'nop = ',(noptotal)
         	write(6,*)'nn = ',nn
         	call exit(-1)
      	endif

!     trick for reading config.in files that are set up as sx,sy,sz and
!     vx,vy,vz; need them in sv(6,nop) matrix
!     it really doesn't have anything to do with old_r arrays, they are
!     just of the appropriate length for this trick

	read(parabundles,*) old_rx,old_ry,old_rz
	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k=1,(noptotal)
         	sv(1,k)=old_rx(k)
         	sv(2,k)=old_ry(k)
         	sv(3,k)=old_rz(k)
      	enddo

      	read(parabundles,*) old_rx,old_ry,old_rz

	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k=1,(noptotal)
         	sv(4,k)=old_rx(k)
         	sv(5,k)=old_ry(k)
         	sv(6,k)=old_rz(k)
      	enddo
      
#else
      	write(6,*)'starting from a previous run'
      	!LR: I changed the box size to 200 for my chosen concentration. This is only valid for 10mM A16
      	boxl = boxlength
      	!open(unit=7,file='results/run'//fname_digits_pre//'.config',status='unknown',form='unformatted') 
 
      	do while (.true.)
		read(precf,end=120) coll,t,old_rx,old_ry,old_rz
      	enddo

120	continue

      !LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k=1,(noptotal)
         	sv(1,k)=old_rx(k)-dnint(old_rx(k))
         	sv(2,k)=old_ry(k)-dnint(old_ry(k))
         	sv(3,k)=old_rz(k)-dnint(old_rz(k))
      	enddo
       	
      	!open(unit=7,file='results/run'//fname_digits_pre//'.lastvel',status='unknown',form='unformatted')
      	read(prelasvel) coll,old_rx,old_ry,old_rz
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k=1,(noptotal)   
         	sv(4,k)=old_rx(k)
         	sv(5,k)=old_ry(k)
         	sv(6,k)=old_rz(k)
      	enddo

#endif

      	!open(unit=7,file='parameters/firstside1.data',status='unknown')
      	read(parafs1,*) fside1
      	l=1

      	do k=1,chnln1 
         	if (fside1(k) .ne. 0) then
	    		fside1(k)=3*chnln1+l
	    		write(6,*) k, fside1(k)
	    		l=l+1
	 	endif
      	enddo

      	write(6,*)' '

	if (chaptype .eq. 1) then
      		!open(unit=7,file='parameters/firstside2.data',status='unknown')
      		read(parafs2,*) fside2
      		!close(unit=7)
      		l=1
      		do k=1,chnln2 
         		if (fside2(k) .ne. 0) then
	    			fside2(k)=nop1+3*chnln2+l
	    			write(6,*) k, fside2(k)
	    			l=l+1
	 		endif
      		enddo
	endif

      	write(6,*)' '
      	!open(unit=7,file='parameters/rn.data',status='unknown')

      	if (chaptype .eq. 1) then
      		do k=1,chnln1+chnln2
         		read(pararn,*) bl_rn(k)
      		enddo
      	else
      		do k=1,chnln1
         		read(pararn,*) bl_rn(k)
      		enddo
      	endif

      	!close(unit=7)
      	!open(unit=7,file='parameters/rc.data',status='unknown')

      	if (chaptype .eq. 1) then
      		do k=1,chnln1+chnln2
         		read(pararc,*) bl_rc(k)
      		enddo
      	else
      		do k=1,chnln1
         		read(pararc,*) bl_rc(k)
      		enddo
      	endif

      	!close(unit=7)
!     	read in sigma, welldia, and epsilon arrays
!     	read in nitrogen bead diameter
      	read(parapro,*)sigma(1)
!     	read in alpha carbon bead diameter 
      	read(parapro,*)sigma(2)
!     	read in side chain bead diameter
      	read(parapro,*)sigma(3)
!     	read in carbonyl carbon bead diameter
      	read(parapro,*)sigma(4)
!     	read in nitrogen well diameter
      	read(parapro,*)welldia(1)
!     	read in dummy well diameter
      	read(parapro,*)welldia(2)
!     	read in side chain well diameter
      	read(parapro,*)welldia(3)
!     	read in carbonyl carbon well diameter
      	read(parapro,*)welldia(4)
!     	read in nitrogen well depth
      	read(parapro,*)epsilon(1)
!     	read in dummy well depth
      	read(parapro,*)epsilon(2)
!     	read in side chain well depth
     	read(parapro,*)epsilon(3)
!     	read in carbonyl carbon well depth
      	read(parapro,*)epsilon(4)

!     	5 through 8 are if bonded, 1 through 4 if not bonded
      	sigma(5)=sigma(1)
      	sigma(6)=sigma(2)
      	sigma(7)=sigma(3)
      	sigma(8)=sigma(4)
      	epsilon(5)=epsilon(1)
      	epsilon(6)=epsilon(2)
      	epsilon(7)=epsilon(3)
      	epsilon(8)=epsilon(4)
      	welldia(5)=welldia(1)
      	welldia(6)=welldia(2)
      	welldia(7)=welldia(3)
      	welldia(8)=welldia(4)
!     	residues in order: grndqehkpstacilmfwyv

      	do k=9,28
		sigma(k)=0.0d0
         	welldia(k)=0.0d0
         	epsilon(k)=0.0d0
      	end do

!      	open(unit=7,file='parameters/identity.inp',status='unknown')
!     	really not reading aa array, just using it temporarily b/c
!     	it's the right length; will re-write it below
!     	input identity is for single chain, need identity(nop)
      	read (paraid,*) aa
!      	close(unit=7)

		write(6,*)'!!!!!!!!!!!aa: ',aa		
		
      	do l=1,nop1,numbeads1
		do k=1,numbeads1
			identity(l+k-1)=aa(k)
		enddo
	enddo

      	do l=nop1+1,nop1+nop2,numbeads2
         	do k=1,numbeads2
            		identity(l+k-1)=aa(numbeads1+k)
         	enddo
      	enddo
      	
      	!LR: This loop sets the identity of the nanoparticle "side-chain"
      	!It is functionally identical to the above loops for the 1st and 2nd species.
!      	open(unit=7,file='parameters/bdln.inp',status='unknown')
      	read(parabdln,*) bdln
!      	close(unit=7)
!      	open(unit=7,file='parameters/hp1.inp',status='unknown')
      	read (parahp1,*) hp1
!      	close(unit=7)
!	open(unit=7,file='parameters/hp2.inp',status='unknown')
      	read (parahp2,*) hp2
!      	close(unit=7)
      	!LR: Reads in the 3rd species parameters
      	

      	do i=1,numbeads1
        	hp(i)=hp1(i)
      	enddo

      	do i=1,numbeads2
       	hp(numbeads1+i)=hp2(i)
      	enddo
      	
      	!LR: Sets the third species parameters

      	do k = 1, nop1
	 	chnnum(k)=(k-1)/numbeads1+1
      	end do

      	do k = 1, nop2
	 	chnnum(nop1+k)=(nop1/numbeads1)+(k-1)/numbeads2+1
      	end do
      	
      	!LR: Sets the "chain number" for the 3rd species
      	!LR: This is a copy of the above loop for the 2nd species, with corrected
      	!range and mathematical conversions
      	!do k = 1, nop3
	 	!chnnum(nop1+nop2+k)=(nop1/numbeads1)+(nop2/numbeads2)+(k-1)/numbeads3+1
      	!end do

#ifdef no_h

	do k=3*chnln1+1,numbeads1
         	hp(k)=0
      	enddo

      	do k=numbeads1+3*chnln2+1,(numbeads1+numbeads2)
         	hp(k)=0
      	enddo

#endif

!      	open(unit=7,file='parametersep/ep19p_ha55a_weakhp.data',status='unknown')
      	write(6,*)'read parametersep/ep19p_ha55a_weakhp.data'

      	do i=1,400
      		read(ep19ha55weak,"(2(2x,i2,2x),f8.3)") iiii,jjjj,eptemp
      		ep(iiii,jjjj)=-eptemp
      	end do
	
!	close(unit=7)
!      	open(unit=7,file='parameters/rcarnrco.data',status='unknown')
         write(6,*)'read parameters/rcarnrco.data'

      	do i=1,20
      		read(simrtoall,"(6(f6.3,2x))") drca(i),drnh(i),drco(i),del_rca(i),del_rnh(i),del_rco(i)
      		if(del_rca(i).lt.del) del_rca(i)=del
      		if(del_rnh(i).lt.del) del_rnh(i)=del
      		if(del_rco(i).lt.del) del_rco(i)=del
      	end do
      
!	close(unit=7)

      	if (chaptype .eq. 1) then
      		do i=1,chnln1
        		if(fside1(i).ne.0) then
        			iii=aa(fside1(i))-8
        			bdln(i)=drca(iii)
        			bl_rn(i)=drnh(iii)
        			bl_rc(i)=drco(iii)
        			del_bdln(i)=del_rca(iii)
        			del_blrn(i)=del_rnh(iii)
        			del_blrc(i)=del_rco(iii)
        		else
        			bdln(i)=drca(1)
        			bl_rn(i)=drnh(1)
        			bl_rc(i)=drco(1)
        			del_bdln(i)=del_rca(1)
        			del_blrn(i)=del_rnh(1)
        			del_blrc(i)=del_rco(1)
        		endif
      		enddo
      		do i=1,chnln2
        		if(fside2(i).ne.0) then
        			iii=identity(fside2(i))-8
        			bdln(chnln1+i)=drca(iii)
        			bl_rn(chnln1+i)=drnh(iii)
        			bl_rc(chnln1+i)=drco(iii)
        			del_bdln(chnln1+i)=del_rca(iii)
        			del_blrn(chnln1+i)=del_rnh(iii)
        			del_blrc(chnln1+i)=del_rco(iii)
        		else
        			bdln(chnln1+i)=drca(1)
        			bl_rn(chnln1+i)=drnh(1)
        			bl_rc(chnln1+i)=drco(1)
        			del_bdln(chnln1+i)=del_rca(1)
        			del_blrn(chnln1+i)=del_rnh(1)
        			del_blrc(chnln1+i)=del_rco(1)
        		endif
      		enddo
      	else
      		do i=1,chnln1
        		if(fside1(i).ne.0) then
        			iii=aa(fside1(i))-8
        			bdln(i)=drca(iii)
        			bl_rn(i)=drnh(iii)
        			bl_rc(i)=drco(iii)
        			del_bdln(i)=del_rca(iii)
        			del_blrn(i)=del_rnh(iii)
        			del_blrc(i)=del_rco(iii)
        		else
        			bdln(i)=drca(1)
        			bl_rn(i)=drnh(1)
        			bl_rc(i)=drco(1)
        			del_bdln(i)=del_rca(1)
        			del_blrn(i)=del_rnh(1)
        			del_blrc(i)=del_rco(1)
        		endif
      		enddo
	endif

!      	open(unit=7,file='parameters/beadwell_ha55a.data',status='unknown')
          write(6,*)'read parameters/beadwell_ha55a.data'

      	do i=1,400
      		read(simwellha55a,711) iiii,jjjj,bdtemp,wltemp
711  		format(2(2x,i2,2x),2(f6.3,2x))
       	bds(iiii,jjjj)=bdtemp
       	wel(iiii,jjjj)=wltemp
      	end do

!	close(unit=7)

	do i=9,28
		sigma(i) = 1.00d0*bds(i,i)
        	welldia(i) = 1.5d0*sigma(i)
      	end do

	sigma(9) = 1.00d0*bds(9,9)-1.2d0
!	open(unit=7,file='parameters/sqz6to10.data',status='old')
         write(6,*)'read parameters/sqz6to10.data'

	do i=1,20
		read(simsqz,'(5(f6.3,2x))') sz8,sz6,sz7,sz9,sz10
		sqz610(1,i+8) = sz6*2/(sigma(i+8)+sigma(4))
       		sqz610(2,i+8) = sz7*2/(sigma(i+8)+sigma(1))
       		sqz610(3,i+8) = sz8*2/(sigma(i+8)+sigma(2))
       		sqz610(4,i+8) = sz9*2/(sigma(i+8)+sigma(2))
       		sqz610(5,i+8) = sz10*2/(sigma(i+8)+sigma(4))
      	end do

!      	close(unit=7)
!      	open(unit=7,file='parameters/mass.data',status='unknown')
          write(6,*)'read parameters/mass.data'

	do i=1,23
       	read(simmass,'(4x,i2,2x,f8.3)') iiii,bmass_temp
		bmass(iiii)=bmass_temp
       end do 

	bmass(3)=bmass(20)

	do i=1,4
		bmass(i+4)=bmass(i)
	end do

!	close(unit=7)

	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do i=1,(noptotal)
		bm(i)=bmass(identity(i))
	end do

!	open(unit=73,file='results/run'//fname_digits_pre//'.rca',status='unknown')
         write(6,*)'read results/.rca file'

	do l=1,nop1,numbeads1
		do k=1,chnln1
			if(fside1(k).ne.0) then
				rx2=(sv(1,l-1+fside1(k))-sv(1,l+k-1))
            			rx2=(rx2-dnint(rx2))**2
            			ry2=(sv(2,l-1+fside1(k))-sv(2,l+k-1))
            			ry2=(ry2-dnint(ry2))**2
            			rz2=(sv(3,l-1+fside1(k))-sv(3,l+k-1))
            			rz2=(rz2-dnint(rz2))**2
            			rca2=dsqrt(rx2+ry2+rz2)
            			write(prerca,7373) "rca ",l,k,aa(k),rca2*boxl,bdln(k) 
            			delx=del_bdln(k)*bdln(k)
            			if(rca2*boxl.lt.(bdln(k)-delx).or.rca2*boxl.gt.(bdln(k)+delx)) then
            				write(prerca,*) "overlap"
            			endif
            			rx2=(sv(1,l-1+fside1(k))-sv(1,l+k+chnln1-1))
            			rx2=(rx2-dnint(rx2))**2
            			ry2=(sv(2,l-1+fside1(k))-sv(2,l+k+chnln1-1))
            			ry2=(ry2-dnint(ry2))**2
            			rz2=(sv(3,l-1+fside1(k))-sv(3,l+k+chnln1-1))
            			rz2=(rz2-dnint(rz2))**2
            			rnh2=dsqrt(rx2+ry2+rz2)
            			write(prerca,7373) "rnh ",l,k,aa(k),rnh2*boxl,bl_rn(k) 
            			delx=del_blrn(k)*bl_rn(k)
            			if(rnh2*boxl.lt.(bl_rn(k)-delx).or.rnh2*boxl.gt.(bl_rn(k)+delx)) then
            				write(prerca,*) "overlap"
            			endif
            			rx2=(sv(1,l-1+fside1(k))-sv(1,l+k+2*chnln1-1))
            			rx2=(rx2-dnint(rx2))**2
            			ry2=(sv(2,l-1+fside1(k))-sv(2,l+k+2*chnln1-1))
            			ry2=(ry2-dnint(ry2))**2
            			rz2=(sv(3,l-1+fside1(k))-sv(3,l+k+2*chnln1-1))
            			rz2=(rz2-dnint(rz2))**2
            			rco2=dsqrt(rx2+ry2+rz2)
            			write(prerca,7373) "rco ",l,k,aa(k),rco2*boxl,bl_rc(k) 
            			delx=del_blrc(k)*bl_rc(k)
            			if(rco2*boxl.lt.(bl_rc(k)-delx).or.rco2*boxl.gt.(bl_rc(k)+delx)) then
            				write(prerca,*) "overlap"
            			endif
7373   			format(a4,2(i4,1x),i2,1x,2(f7.4,1x))
        		endif
		end do
	end do

	if (chaptype .eq. 1) then
      		do l=1,nop2,numbeads2
        		do k=1,chnln2
        			if(fside2(k).ne.0) then
            				rx2=(sv(1,l-1+fside2(k))-sv(1,nop1+l+k-1))
            				rx2=(rx2-dnint(rx2))**2
            				ry2=(sv(2,l-1+fside2(k))-sv(2,nop1+l+k-1))
            				ry2=(ry2-dnint(ry2))**2
            				rz2=(sv(3,l-1+fside2(k))-sv(3,nop1+l+k-1))
            				rz2=(rz2-dnint(rz2))**2
            				rca2=dsqrt(rx2+ry2+rz2)
            				write(prerca,7373) "rca ",l,k,aa(k),rca2*boxl,bdln(chnln1+k) 
            				delx=del_bdln(chnln1+k)*bdln(chnln1+k)
            				if(rca2*boxl.lt.(bdln(chnln1+k)-delx).or.rca2*boxl.gt.(bdln(chnln1+k)+delx)) then
            					write(prerca,*) "overlap"
            				endif
            				rx2=(sv(1,l-1+fside2(k))-sv(1,nop1+l+k+chnln2-1))
            				rx2=(rx2-dnint(rx2))**2
            				ry2=(sv(2,l-1+fside2(k))-sv(2,nop1+l+k+chnln2-1))
            				ry2=(ry2-dnint(ry2))**2
            				rz2=(sv(3,l-1+fside2(k))-sv(3,nop1+l+k+chnln2-1))
            				rz2=(rz2-dnint(rz2))**2
            				rnh2=dsqrt(rx2+ry2+rz2)
            				write(prerca,7373) "rnh ",l,k,aa(k),rnh2*boxl,bl_rn(chnln1+k) 
            				delx=del_blrn(chnln1+k)*bl_rn(chnln1+k)
            				if(rnh2*boxl.lt.(bl_rn(chnln1+k)-delx).or.rnh2*boxl.gt.(bl_rn(chnln1+k)+delx)) then
            					write(prerca,*) "overlap"
            				endif
            				rx2=(sv(1,l-1+fside2(k))-sv(1,nop1+l+k+2*chnln2-1))
            				rx2=(rx2-dnint(rx2))**2
            				ry2=(sv(2,l-1+fside2(k))-sv(2,nop1+l+k+2*chnln2-1))
            				ry2=(ry2-dnint(ry2))**2
            				rz2=(sv(3,l-1+fside2(k))-sv(3,nop1+l+k+2*chnln2-1))
            				rz2=(rz2-dnint(rz2))**2
            				rco2=dsqrt(rx2+ry2+rz2)
            				write(prerca,7373) "rco ",l,k,aa(k),rco2*boxl,bl_rc(chnln1+k) 
            				delx=del_blrc(chnln1+k)*bl_rc(chnln1+k)
            				if(rco2*boxl.lt.(bl_rc(chnln1+k)-delx).or.rco2*boxl.gt.(bl_rc(chnln1+k)+delx)) then
            					write(prerca,*) "overlap"
            				endif
				endif
        		end do
      		end do
      	endif

!      	close(unit=73)
!     	scale velocities to desired temperature
      	sumvel=0.d0

	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
     	do k=1,(noptotal)
		sumvel=sumvel+bm(k)*(sv(4,k)*sv(4,k)+sv(5,k)*sv(5,k)+sv(6,k)*sv(6,k))
      	enddo

		!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	tred=sumvel/3.d0/dble(noptotal)
      	write(6,*)'initial temperature = ',tred
      	write(6,*)'desired temperature = ',setemp

#ifdef ctemp

      	sumx=0.0
      	sumy=0.0
      	sumz=0.0  
       
       !LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k=1,(noptotal)         
         	x1=drandm(0)*0.75d0+0.2d0
		x2=drandm(0)*0.75d0+0.2d0
         	y1=drandm(0)*0.75d0+0.2d0
         	y2=drandm(0)*0.75d0+0.2d0
         	sv(4,k)=dsqrt(-2.d0*log(x1))*cos(2.d0*pi*y1)
         	sv(5,k)=dsqrt(-2.d0*log(x1))*sin(2.d0*pi*y1)
         	sv(6,k)=dsqrt(-2.d0*log(x2))*cos(2.d0*pi*y2)
      	end do
      
!     	scale velocities so that linear momentum = 0

		!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k=1,(noptotal)
         	sumx=sumx+sv(4,k)
         	sumy=sumy+sv(5,k)
         	sumz=sumz+sv(6,k)
      	end do
      
      !LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k=1,(noptotal)         
         	sv(4,k)=sv(4,k)-sumx/real((noptotal))
         	sv(5,k)=sv(5,k)-sumy/real((noptotal))
         	sv(6,k)=sv(6,k)-sumz/real((noptotal))
      	end do    
 
      	sumvel=0.d0

		!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k=1,(noptotal)
		sumvel=sumvel+(sv(4,k)*sv(4,k)+sv(5,k)*sv(5,k)+sv(6,k)*sv(6,k))/bm(k)
      	enddo  ! since here sv is actually momentum not velocity

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	tred=sumvel/3.d0/dble((noptotal))
      	const=dsqrt(setemp/tred)

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k=1,(noptotal)
		sv(4,k)=const*sv(4,k)/bm(k)
         	sv(5,k)=const*sv(5,k)/bm(k)
         	sv(6,k)=const*sv(6,k)/bm(k)
      	enddo

      	sumvel=0.d0

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	do k=1,(noptotal)
         	sumvel=sumvel+bm(k)*(sv(4,k)*sv(4,k)+sv(5,k)*sv(5,k)+sv(6,k)*sv(6,k))
      	enddo

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
      	tred=sumvel/3.d0/dble((noptotal))
      	write(6,*)'new scaled temperature = ',tred

#endif

      	write(6,*)
      	write(6,*)'***********************************************'
     	write(6,*)
     	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
     	write(6,*)'number of particles = ',(noptotal)
     	write(6,*)'set temperature = ',setemp
      	write(6,*)'time interval of false steps = ', interval
      	write(6,*)'chnln1 and chnln2 = ',chnln1,chnln2 
      	write(6,*)'numbeads1 and numbeads1 = ',numbeads1,numbeads2
      	write(6,*)
      	write(6,*)'n_b_hbond=',n_b_hbond
      	write(6,*)'n_b_hydro=',n_b_hydro
      	write(6,*)
      	write(6,*)
      	write(6,*)'nitrogen bead diameter = ',sigma(1)
      	write(6,*)'alpha carbon bead diameter = ',sigma(2)
      	write(6,*)'carbonyl carbon bead diameter = ',sigma(4) 
      	write(6,*)
      	write(6,*)'nitrogen well diameter = ',welldia(1)
      	write(6,*)'alpha carbon well  diameter = ',welldia(2)
      	write(6,*)'carbonyl carbon well  diameter = ',welldia(4)
      	write(6,*)
      	write(6,*)'nitrogen well depth = ',epsilon(1)
      	write(6,*)'alpha carbon well depth = ',epsilon(2)
      	write(6,*)'side chain well depth = ',epsilon(3)
      	write(6,*)'carbonyl carbon well depth = ',epsilon(4)
      	write(6,*)
      	write(6,*)'n - c alpha bond length = ',dnc
      	write(6,*)'c alpha - c bond length = ',dcc
      	write(6,*)'c alpha - n bond length = ',dcn
      	write(6,*)'c alpha - side chain bond length = (variable)'
      	write(6,*)
88   	format(i5,2f8.5)
      	write(6,*)'amino acid sequence is:'
      	write(6,*)aa
      	write(6,*) 
      	write(6,*)'***********************************************'
      	write(6,*)   
      	write(6,*)'sc idents, sizes,? well widths for chain 1:'

      	do k=3*chnln1+1,numbeads1
		write(6,88)identity(k),sigma(identity(k)), &
&		welldia(identity(k))
      	enddo

      	write(6,*)
      	write(6,*)'sc idents, sizes,? well widths for chain 2:'

      	do k=nop1+3*chnln2+1,nop1+numbeads2
		write(6,88)identity(k),sigma(identity(k)), &
&		welldia(identity(k))
      	enddo

      	write(6,*)
      	write(6,*)'***********************************************'
      	write(6,*)

      	do k=1,numbeads1
		write(6,89) k,identity(k),bm(k)   ! by mookyung
      	end do

      	do k=nop1+1,nop1+numbeads2
		write(6,89) k,identity(k),bm(k)   ! by mookyung
      	end do
      	
      	!LR: Prints the identity of the 3rd species bead

89	format(i6,2x,i5,2f8.5)
      	write(6,*)
      	write(6,*)'***********************************************'
      	write(6,*)   
      	write(6,*)'hp sequence is:'
      	write(6,*) hp
      	write(6,*)   
      	write(6,*)'***********************************************'

      	do i=9,28
      		write(6,'(5(f6.3,2x))') (sqz610(k,i),k=1,5)
      	end do

      	write(6,*)   

      	call scale_down()

      	return

      	end




