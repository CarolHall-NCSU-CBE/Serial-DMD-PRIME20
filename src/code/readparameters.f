subroutine readparameters
	use global
	use inputreadin
	implicit none
	
	integer aa(numbeads1+numbeads2),i,iiii,jjjj,l,k
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

    	
end subroutine
