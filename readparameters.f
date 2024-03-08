subroutine readparameters
	use global
	use inputreadin
	implicit none
	
	integer aa(numbeads1+numbeads2),i,iiii,jjjj,l,k,kk
	real*8 eptemp,bdtemp,wltemp,bmass_temp !! by mookyung
	
	call chdir(rundir)
	open(unit=paraid,file='parameters/identity.inp',status='unknown')
	read (paraid,*) aa
	close(paraid)	
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
	open(unit=parahp1,file='parameters/hp1.inp',status='unknown')
      	read (parahp1,*) hp1
      	close(parahp1)
	open(parahp2,file='parameters/hp2.inp',status='unknown')
      	read (parahp2,*) hp2
	close(parahp2)
      	do i=1,numbeads1
        	hp(i)=hp1(i)
      	enddo
      	do i=1,numbeads2
       		hp(numbeads1+i)=hp2(i)
      	enddo
	open(unit=parafs1,file='parameters/firstside1.data',status='unknown')
	read(parafs1,*) fside1
      	close(parafs1)
	l=1
      	do k=1,chnln1 
        	if (fside1(k) .ne. 0) then
	    		fside1(k)=3*chnln1+l
	    		l=l+1
	 	endif
      	enddo	
	if (chaptype .eq. 1) then
		open(unit=parafs2,file='parameters/firstside2.data',status='unknown')
      	      	read(parafs2,*) fside2
      		close(parafs2)
      		l=1
      		do k=1,chnln2 
         		if (fside2(k) .ne. 0) then
	    			fside2(k)=nop1+3*chnln2+l
	    			l=l+1
	 		endif
      		enddo
	endif
	call chdir(mydir) 
	open(unit=parapro,file='parameters/protein.data',status='unknown')
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
	close(parapro)
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
	open(unit=simwellha55a,file='parameters/beadwell_ha55a.data',status='unknown')
      	do i=1,400
      		read(simwellha55a,711) iiii,jjjj,bdtemp,wltemp
711  		format(2(2x,i2,2x),2(f6.3,2x))
       		bds(iiii,jjjj)=bdtemp
       		wel(iiii,jjjj)=wltemp
      	end do
	close(simwellha55a)
	do k = 9,28
         	do kk = 9,28
            		welldia_sq(k,kk) = wel(k,kk)**2/boxl**2
         	enddo
   	enddo

!     	do k=9,28
!		sigma(i) = 1.00d0*bds(i,i)
!       	welldia(i) = 1.5d0*sigma(i)
!       	epsilon(k)=0.0d0
!    	end do
    	
end subroutine