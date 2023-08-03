subroutine genconfig
	use inputreadin
	use global
     	implicit none
        integer nb1p,nb21p,chnln1p

      parameter(nb1p=124)       ! number of beads 124, we don't change lines 6-8
      parameter(nb21p=124)       ! number of beads 124	  
      parameter(chnln1p=31)     ! chain length  31

	integer aa(nb1),identity1(nopwg1),bptnr1(nop1)  
	integer aa2(nb2),identity2(nopwg2),bptnr2(nop2)	  
      integer i,j,k,l,iii,iiii,numbeads,jjjj,ii
	real*8 xr
      real*8 old_rxgenpep(nb1p), old_rygenpep(nb1p), old_rzgenpep(nb1p)     
      real*8 dave_rx(124), dave_ry(124), dave_rz(124)     	  
      real*8 new_rx(numbeads1), new_ry(numbeads1), new_rz(numbeads1)
      real*8 new_rx2(numbeads2), new_ry2(numbeads2), new_rz2(numbeads2)	  
      real*8 new_rxf(numbeads1), new_ryf(numbeads1), new_rzf(numbeads1) 	  
      real*8 new_rxf2(numbeads2), new_ryf2(numbeads2), new_rzf2(numbeads2)
      real*8 xtmp(nop1+nop2),ytmp(nop1+nop2),ztmp(nop1+nop2)            
      real*8 svgcf(6,nopwg1+nopwg2),boxlorig
      real*8 xmin,xmax,ymin,ymax,zmin,zmax
      real*8 dislocx(numbeads1),dislocy(numbeads1),dislocz(numbeads1)
      real*8 dislocx2(numbeads2),dislocy2(numbeads2),dislocz2(numbeads2)
      real sig_ij,sig_max,bdsgencf(400,400),welcf(400,400),bdtemp,wltemp
 
      real*8 rmin,r3,rca,rnh,rco,DS1, DS2, DS3                              
      real*8 rcax,rcay,rcaz,rnhx,rnhy,rnhz,rcox,rcoy,rcoz
      real*8 reducedtemp !del, pi, setemp
      real*8 drca(20),drnh(20),drco(20)
      real*8 del_rca(20),del_rnh(20),del_rco(20)
      real*8 bdln1(chnln1),bl_rn1(chnln1),bl_rc1(chnln1),del_bdln1(chnln1),del_blrn1(chnln1),del_blrc1(chnln1)
      real*8 bdln2(chnln2),bl_rn2(chnln2),bl_rc2(chnln2),del_bdln2(chnln2),del_blrn2(chnln2),del_blrc2(chnln2)
      real*8 scalerca,rx2,ry2,rz2
      real*8 d1,d2,d3,d1min,d2min,d3min,move(chnln1),move2(chnln2)
      real*8 rxij,ryij,rzij,rijsq
      REAL*8 sqz610cf(5,20)
      REAL*8 SZ6,SZ7,SZ8,SZ9,SZ10
      REAL*8 SZ6X,SZ7X,SZ8X,SZ9X,SZ10X
      REAL*8 SZ6Y,SZ7Y,SZ8Y,SZ9Y,SZ10Y
      REAL*8 SZ6Z,SZ7Z,SZ8Z,SZ9Z,SZ10Z
      INTEGER ISZ,ID
	  
      real*8 bm1(nopwg1),bm2(nopwg2),bmass_temp,bmass(28)                     
      real*8 sumvel,tred,const,sumx,sumy,sumz,x1,x2,y1,y2
	real redsimT  

	integer peplength1, peplength2, resid
	character rescf, residue(28), sc

	!character*517 ::  path, mydir, rundir
	!integer realpath
	!logical :: back=.true.
	
	logical exist_flag1,exist_flag2
      inquire(file = 'inputs/chninfo-n1.data', exist = exist_flag1)	 
      inquire(file = 'inputs/chninfo-n2.data', exist = exist_flag2)	  
      iflag = 1331171207 !1058472402   !197071101   !1331171207
      xr=drandm(iflag)
      call srand(iflag)
	boxl = boxlength
	reducedtemp=0.50
      setemp = reducedtemp*12

	call gencf_file_opener() 
      read(genericpepx,*) numbeads,boxlorig,old_rxgenpep
      read(genericpepy,*) numbeads,boxlorig,old_rygenpep  
      read(genericpepz,*) numbeads,boxlorig,old_rzgenpep

7       format(A4,3X,I4,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3)
	  
!      open(7,file='ala31.pdb',status='unknown',form='formatted')
!      do i=1,31
!      write(7,7) 'ATOM',(i-1)*4+1,' N  ','ALA','A',i,dave_rx(i+31),dave_ry(i+31),dave_rz(i+31)
!      write(7,7) 'ATOM',(i-1)*4+2,' CA ','ALA','A',i,dave_rx(i),dave_ry(i),dave_rz(i)
!      write(7,7) 'ATOM',(i-1)*4+3,' C  ','ALA','A',i,dave_rx(i+2*31),dave_ry(i+2*31),dave_rz(i+2*31)
!      write(7,7) 'ATOM',(i-1)*4+4,' CB ','ALA','A',i,dave_rx(i+3*31),dave_ry(i+3*31),dave_rz(i+3*31)
!      end do
!      close(7)	  
	  
!      open(7,file='abeta1-40.pdb',status='old')
!      do i=1,chnln1
!      read(7,'(30x,3(f8.3))') old_rxgenpep(i+chnln1),old_rygenpep(i+chnln1),old_rzgenpep(i+chnln1)
!      read(7,'(30x,3(f8.3))') old_rxgenpep(i),old_rygenpep(i),old_rzgenpep(i)
!      read(7,'(30x,3(f8.3))') old_rxgenpep(i+2*chnln1),old_rygenpep(i+2*chnln1),old_rzgenpep(i+2*chnln1)
!      read(7,'(30x,3(f8.3))') old_rxgenpep(i+3*chnln1),old_rygenpep(i+3*chnln1),old_rzgenpep(i+3*chnln1)
!      end do
!      close(7)
 
! Create Input file - VN

	do i=1,23
       		read(massfile,*) rescf, resid, bmass_temp
		if (resid .ge. 9) then
			residue(resid:resid)=rescf
		endif
		bmass(resid) = bmass_temp
       	end do
	peplength1 = len_trim(pep1)
	peplength2 = len_trim(pep2)

	do i = 1,peplength1
		write(ideninp,*) 2
	enddo
	do i = peplength1+1, 2*peplength1
		write(ideninp,*) 1
	enddo
	do i = 2*peplength1+1,3*peplength1
		write(ideninp,*) 4
	enddo
	do i = 3*peplength1+1, 4*peplength1
		write(ideninp,*) maxloc(scan(residue,pep1((i-3*peplength1):(i-3*peplength1))),1)
	enddo

	do i = 1,peplength2
		write(iden2inp,*) 2
	enddo
	do i = peplength2+1, 2*peplength2
		write(iden2inp,*) 1
	enddo
	do i = 2*peplength2+1,3*peplength2
		write(iden2inp,*) 4
	enddo
	do i = 3*peplength2+1, 4*peplength2
		write(iden2inp,*) maxloc(scan(residue,pep2((i-3*peplength2):(i-3*peplength2))),1)
	enddo

!read in identities
      rewind(ideninp)
      read(ideninp,*) aa
      do l=1,nopwg1,nb1
       do k=1,nb1 
		identity1(l+k-1)=aa(k)
      enddo
      enddo
	
	rewind(iden2inp)
      read(iden2inp,*) aa2    
      do l=1,nopwg2,nb2
       do k=1,nb2 
		identity2(l+k-1)=aa2(k)
      enddo
      enddo	  
	  
	!write file to verify identities
      !open(checkgencf,file='checks/identity.out')
      do i=1,nopwg1
		write(checkid,*) i,identity1(i)
      enddo
      do i=1,nopwg2
		write(checkid,*) i,identity2(i)
      enddo	  
      !close(7)
	
	!this writes the hydrophobicity inputs
      !open(7,file='parameters/hp1.inp')
      do i=1,nb1
      if (identity1(i) .le. 4) then
			write(parahp1,*)'0'
      elseif (identity1(i) .gt. 9) then
			write(parahp1,*)'1'			
      endif
      enddo
      !close(7)
      !open(7,file='parameters/hp2.inp')
      do i=1,nb2
      if (identity2(i) .le. 4) then
		write(parahp2,*)'0'
      elseif (identity2(i) .gt. 9) then
		write(parahp2,*)'1'			
      endif
      enddo
      !close(7)

	!this writes the sidechain inputs
      !open(7,file='parameters/firstside1.data')
      do i=chnln1*3+1,nb1
      if (identity1(i) .eq. 9) then
		write(parafs1,*)'0'
      else
		write(parafs1,*)'1'
      endif
      enddo
      !close(7)
      !open(7,file='parameters/firstside2.data')
      do i=chnln2*3+1,nb2
      if (identity2(i) .eq. 9) then
		write(parafs2,*)'0'
      else
		write(parafs2,*)'1'
      endif
      enddo
      !close(7)
	  
	!this writes the identity input for the simulation
      !open(7,file='parameters/identity.inp')   
      do l=1,nb1
      	if (identity1(l) .ne. 9) write(paraid,*) identity1(l)
      enddo
      do l=1,nb2
      	if (identity2(l) .ne. 9) write(paraid,*) identity2(l)
      enddo
      !close(7)

	!del = 0.02375d0

	!read in l/d-isomer pseudobond constraints
	!call chdir(mydir)
      !open(7,file='genconfig/parameters/rcarnrco.data',status='old')
      do i=1,20
		read(parartoall,"(6(f6.3,2x))") drca(i),drnh(i),drco(i),del_rca(i),del_rnh(i),del_rco(i)
	      if(del_rca(i).lt.del) del_rca(i)=del
      		if(del_rnh(i).lt.del) del_rnh(i)=del
      		if(del_rco(i).lt.del) del_rco(i)=del
      end do
      !close(7)

	!reassign pseudobond identities
      do i=1,chnln1
		iii=identity1(3*chnln1+i)-8
        	bdln1(i)=drca(iii)
        	bl_rn1(i)=drnh(iii)
        	bl_rc1(i)=drco(iii)
        	del_bdln1(i)=del_rca(iii)
        	del_blrn1(i)=del_rnh(iii)
        	del_blrc1(i)=del_rco(iii)
      end do

      do i=1,chnln2
		iii=identity2(3*chnln2+i)-8
        	bdln2(i)=drca(iii)
        	bl_rn2(i)=drnh(iii)
        	bl_rc2(i)=drco(iii)
        	del_bdln2(i)=del_rca(iii)
        	del_blrn2(i)=del_rnh(iii)
        	del_blrc2(i)=del_rco(iii)
      end do

       !OPEN(UNIT=7,file='genconfig/parameters/sqz6to10.data',STATUS='OLD')
      DO I=1,20
       READ(parasqz,'(5(f6.3,2x))') (sqz610cf(K,I),K=1,5)
      END DO
       !CLOSE(UNIT=7)
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      if((exist_flag1).and.(exist_flag2)) goto 999	! Yiming
	
        do i=1,chnln1
        	move(i) = dadjust1 
        enddo
	
         !this entire section takes each sidechain and moves the location of the bead in a box until it satisfies the pseudobond contraints
      do k=3*chnln1p+1,3*chnln1p+chnln1
		rx2=(old_rxgenpep(k)-old_rxgenpep(k-3*chnln1p))
		ry2=(old_rygenpep(k)-old_rygenpep(k-3*chnln1p))
		rz2=(old_rzgenpep(k)-old_rzgenpep(k-3*chnln1p))
		scalerca=bdln1(k-3*chnln1p)/dsqrt(rx2**2+ry2**2+rz2**2)
      		old_rxgenpep(k)=old_rxgenpep(k-3*chnln1p)+scalerca*rx2
      		old_rygenpep(k)=old_rygenpep(k-3*chnln1p)+scalerca*ry2
      		old_rzgenpep(k)=old_rzgenpep(k-3*chnln1p)+scalerca*rz2
      end do
                ds1 = 0.02
                ds2 = 0.02
                ds3 = 0.02
      write(6,*)
      write(6,'(a50)') "first value should be less than the 0.01 tolerance"
      write(6,'(a62)') "if not, adjust d1 d2 d3 values in sidechain adjustment section"

      do k=1,chnln1
		rmin=100.d0		
         	do d1=-move(k),move(k),0.005d0
			rcax=old_rxgenpep(k+3*chnln1p)+d1-old_rxgenpep(k)
          		rnhx=old_rxgenpep(k+3*chnln1p)+d1-old_rxgenpep(k+chnln1p)
          		rcox=old_rxgenpep(k+3*chnln1p)+d1-old_rxgenpep(k+2*chnln1p)
             if(k.ge.2) sz6x=old_rxgenpep(k+3*chnln1p)+d1-old_rxgenpep(k-1+2*chnln1p)
             if(k.ge.2) sz8x=old_rxgenpep(k+3*chnln1p)+d1-old_rxgenpep(k-1)
             if(k.lt.chnln1) sz7x=old_rxgenpep(k+3*chnln1p)+d1-old_rxgenpep(k+1+chnln1p)
             if(k.lt.chnln1) sz9x=old_rxgenpep(k+3*chnln1p)+d1-old_rxgenpep(k+1)
         		do d2=-move(k),move(k),.005d0
           			rcay=old_rygenpep(k+3*chnln1p)+d2-old_rygenpep(k)
           			rnhy=old_rygenpep(k+3*chnln1p)+d2-old_rygenpep(k+chnln1p)
           			rcoy=old_rygenpep(k+3*chnln1p)+d2-old_rygenpep(k+2*chnln1p)
             if(k.ge.2) sz6y=old_rygenpep(k+3*chnln1p)+d2-old_rygenpep(k-1+2*chnln1p)
             if(k.ge.2) sz8y=old_rygenpep(k+3*chnln1p)+d2-old_rygenpep(k-1)
             if(k.lt.chnln1) sz7y=old_rygenpep(k+3*chnln1p)+d2-old_rygenpep(k+1+chnln1p)
             if(k.lt.chnln1) sz9y=old_rygenpep(k+3*chnln1p)+d2-old_rygenpep(k+1)
         			do d3=-move(k),move(k),.005d0
           				rcaz=old_rzgenpep(k+3*chnln1p)+d3-old_rzgenpep(k)
           				rnhz=old_rzgenpep(k+3*chnln1p)+d3-old_rzgenpep(k+chnln1p)
           				rcoz=old_rzgenpep(k+3*chnln1p)+d3-old_rzgenpep(k+2*chnln1p)
             if(k.ge.2) sz6z=old_rzgenpep(k+3*chnln1p)+d3-old_rzgenpep(k-1+2*chnln1p)
             if(k.ge.2) sz8z=old_rzgenpep(k+3*chnln1p)+d3-old_rzgenpep(k-1)
             if(k.lt.chnln1) sz7z=old_rzgenpep(k+3*chnln1p)+d3-old_rzgenpep(k+1+chnln1p)
             if(k.lt.chnln1) sz9z=old_rzgenpep(k+3*chnln1p)+d3-old_rzgenpep(k+1)
           				rca=dsqrt(rcax**2+rcay**2+rcaz**2)
           				rnh=dsqrt(rnhx**2+rnhy**2+rnhz**2)
           				rco=dsqrt(rcox**2+rcoy**2+rcoz**2)
             sz6=dsqrt(sz6x**2+sz6y**2+sz6z**2)
             sz7=dsqrt(sz7x**2+sz7y**2+sz7z**2)
             sz8=dsqrt(sz8x**2+sz8y**2+sz8z**2)
             sz9=dsqrt(sz9x**2+sz9y**2+sz9z**2)
!            sz10=dsqrt(sz10x**2+sz10y**2+sz10z**2)
                ds1 = dabs(bdln1(k)-rca)
                ds2 = dabs(bl_rn1(k)-rnh)
                ds3 = dabs(bl_rc1(k)-rco)
           				r3=dabs(bdln1(k)-rca)+dabs(bl_rn1(k)-rnh)+dabs(bl_rc1(k)-rco)
						 isz=0
             if(k.ge.2) then
             if(sz6.lt.sqz610cf(2,identity1(3*chnln1+k)-8)) isz=1
             if(sz8.lt.sqz610cf(1,identity1(3*chnln1+k)-8)) isz=1
             endif
             if(k.lt.chnln1) then
             if(sz7.lt.sqz610cf(3,identity1(3*chnln1+k)-8)) isz=1
             if(sz9.lt.sqz610cf(4,identity1(3*chnln1+k)-8)) isz=1
             endif
!            if(k.ge.3) then
!            if(sz10.lt.sqz610cf(5,identity1(3*chnln1+k)-8)) isz=1
!            endif
             if(isz.eq.0) id=id+1
             if((r3.lt.rmin).and.(isz.eq.0)) then
            					d1min=d1
          					    d2min=d2
               					d3min=d3
               					rmin=r3
           				endif
         			end do        
         		end do        
      end do        
		
       		old_rxgenpep(k+3*chnln1p)=old_rxgenpep(k+3*chnln1p)+d1min
       		old_rygenpep(k+3*chnln1p)=old_rygenpep(k+3*chnln1p)+d2min
       		old_rzgenpep(k+3*chnln1p)=old_rzgenpep(k+3*chnln1p)+d3min
       		rcax=old_rxgenpep(k+3*chnln1p)-old_rxgenpep(k)
       		rnhx=old_rxgenpep(k+3*chnln1p)-old_rxgenpep(k+chnln1p)
       		rcox=old_rxgenpep(k+3*chnln1p)-old_rxgenpep(k+2*chnln1p)
       		rcay=old_rygenpep(k+3*chnln1p)-old_rygenpep(k)
       		rnhy=old_rygenpep(k+3*chnln1p)-old_rygenpep(k+chnln1p)
       		rcoy=old_rygenpep(k+3*chnln1p)-old_rygenpep(k+2*chnln1p)
       		rcaz=old_rzgenpep(k+3*chnln1p)-old_rzgenpep(k)
       		rnhz=old_rzgenpep(k+3*chnln1p)-old_rzgenpep(k+chnln1p)
       		rcoz=old_rzgenpep(k+3*chnln1p)-old_rzgenpep(k+2*chnln1p)
       		rca=dsqrt(rcax**2+rcay**2+rcaz**2)
       		rnh=dsqrt(rnhx**2+rnhy**2+rnhz**2)
       		rco=dsqrt(rcox**2+rcoy**2+rcoz**2)
             if(k.ge.2) sz6x=old_rxgenpep(k+3*chnln1)-old_rxgenpep(k-1+2*chnln1)
             if(k.ge.2) sz8x=old_rxgenpep(k+3*chnln1)-old_rxgenpep(k-1)
             if(k.lt.chnln1) sz7x=old_rxgenpep(k+3*chnln1)-old_rxgenpep(k+1+chnln1)
             if(k.lt.chnln1) sz9x=old_rxgenpep(k+3*chnln1)-old_rxgenpep(k+1)
             if(k.ge.3) sz10x=old_rxgenpep(k+3*chnln1)-old_rxgenpep(k-2+2*chnln1)
             if(k.ge.2) sz6y=old_rygenpep(k+3*chnln1)-old_rygenpep(k-1+2*chnln1)
             if(k.ge.2) sz8y=old_rygenpep(k+3*chnln1)-old_rygenpep(k-1)
             if(k.lt.chnln1) sz7y=old_rygenpep(k+3*chnln1)-old_rygenpep(k+1+chnln1)
             if(k.lt.chnln1) sz9y=old_rygenpep(k+3*chnln1)-old_rygenpep(k+1)
             if(k.ge.3) sz10y=old_rygenpep(k+3*chnln1)-old_rygenpep(k-2+2*chnln1)
             if(k.ge.2) sz6z=old_rzgenpep(k+3*chnln1)-old_rzgenpep(k-1+2*chnln1)
             if(k.ge.2) sz8z=old_rzgenpep(k+3*chnln1)-old_rzgenpep(k-1)
             if(k.lt.chnln1) sz7z=old_rzgenpep(k+3*chnln1)-old_rzgenpep(k+1+chnln1)
             if(k.lt.chnln1) sz9z=old_rzgenpep(k+3*chnln1)-old_rzgenpep(k+1)
             if(k.ge.3) sz10z=old_rzgenpep(k+3*chnln1)-old_rzgenpep(k-2+2*chnln1)
             sz6=dsqrt(sz6x**2+sz6y**2+sz6z**2)
             sz7=dsqrt(sz7x**2+sz7y**2+sz7z**2)
             sz8=dsqrt(sz8x**2+sz8y**2+sz8z**2)
             sz9=dsqrt(sz9x**2+sz9y**2+sz9z**2)
             sz10=dsqrt(sz10x**2+sz10y**2+sz10z**2)
        	write(6,'(i3,7(f8.4,1x))') k,rmin,d1min,d2min,d3min,rca,rnh,rco
		if (rmin .ge. 0.01) then
			write(6,*)'error in residue',k
			stop 'Configurating generation is terminated'
	        endif
      end do

      pi = 4.d0*datan(1.d0)
      coll = 0
      t=0.d0
  
	!this gets rid of glycine residues
        	l=0
    	 do i=0,3
    		do k=1,chnln1
        		if(identity1(i*chnln1+k).ne.9) then
         			l=l+1
         			new_rx(l)=old_rxgenpep(i*chnln1p+k)
         			new_ry(l)=old_rygenpep(i*chnln1p+k)
         			new_rz(l)=old_rzgenpep(i*chnln1p+k)
              endif
       		enddo
  	    enddo
	!call chdir(rundir)	
      !open(7,file='inputs/chninfo-n1.data',status='unknown',form='formatted')
      do i=1,numbeads1
	      write(inpinfon1,*) new_rx(i),new_ry(i),new_rz(i) ! xyz is stored in 4 groups
      enddo
      !close(unit=7)
	  
	!this makes a pdb to show the structure of one of each species
      !open(checkcf1,file='checks/configone1.pdb')
      write(checkcf1,*) numbeads1
      do i=1,numbeads1
       write(checkcf1,'(a6,i5,a3,1x,15x,3f8.3)') 'atom  ',i,'n', new_rx(i), new_ry(i), new_rz(i)
      enddo
      !close(7)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=1,chnln2
        	move2(i) = dadjust2 
        enddo
	
    	d1min=0.0
        d2min=0.0
     	d3min=0.0
         rmin=0.0
      rewind(genericpepx)  
      read(genericpepx,*) numbeads,boxlorig,old_rxgenpep
	rewind(genericpepy)
      read(genericpepy,*) numbeads,boxlorig,old_rygenpep
	rewind(genericpepz)
      read(genericpepz,*) numbeads,boxlorig,old_rzgenpep
	
         !this entire section takes each sidechain and moves the location of the bead in a box until it satisfies the pseudobond contraints
      do k=3*chnln1p+1,3*chnln1p+chnln2
		rx2=(old_rxgenpep(k)-old_rxgenpep(k-3*chnln1p))
		ry2=(old_rygenpep(k)-old_rygenpep(k-3*chnln1p))
		rz2=(old_rzgenpep(k)-old_rzgenpep(k-3*chnln1p))
		scalerca=bdln2(k-3*chnln1p)/dsqrt(rx2**2+ry2**2+rz2**2)
      		old_rxgenpep(k)=old_rxgenpep(k-3*chnln1p)+scalerca*rx2
      		old_rygenpep(k)=old_rygenpep(k-3*chnln1p)+scalerca*ry2
      		old_rzgenpep(k)=old_rzgenpep(k-3*chnln1p)+scalerca*rz2
      end do
                ds1 = 0.02
                ds2 = 0.02
                ds3 = 0.02
      write(6,*)
      write(6,'(a50)') "first value should be less than the 0.01 tolerance"
      write(6,'(a62)') "if not, adjust d1 d2 d3 values in sidechain adjustment section"

      do k=1,chnln2
		rmin=100.d0		
         	do d1=-move2(k),move2(k),0.005d0
			rcax=old_rxgenpep(k+3*chnln1p)+d1-old_rxgenpep(k)
          		rnhx=old_rxgenpep(k+3*chnln1p)+d1-old_rxgenpep(k+chnln1p)
          		rcox=old_rxgenpep(k+3*chnln1p)+d1-old_rxgenpep(k+2*chnln1p)
             if(k.ge.2) sz6x=old_rxgenpep(k+3*chnln1p)+d1-old_rxgenpep(k-1+2*chnln1p)
             if(k.ge.2) sz8x=old_rxgenpep(k+3*chnln1p)+d1-old_rxgenpep(k-1)
             if(k.lt.chnln2) sz7x=old_rxgenpep(k+3*chnln1p)+d1-old_rxgenpep(k+1+chnln1p)
             if(k.lt.chnln2) sz9x=old_rxgenpep(k+3*chnln1p)+d1-old_rxgenpep(k+1)
         		do d2=-move2(k),move2(k),.005d0
           			rcay=old_rygenpep(k+3*chnln1p)+d2-old_rygenpep(k)
           			rnhy=old_rygenpep(k+3*chnln1p)+d2-old_rygenpep(k+chnln1p)
           			rcoy=old_rygenpep(k+3*chnln1p)+d2-old_rygenpep(k+2*chnln1p)
             if(k.ge.2) sz6y=old_rygenpep(k+3*chnln1p)+d2-old_rygenpep(k-1+2*chnln1p)
             if(k.ge.2) sz8y=old_rygenpep(k+3*chnln1p)+d2-old_rygenpep(k-1)
             if(k.lt.chnln2) sz7y=old_rygenpep(k+3*chnln1p)+d2-old_rygenpep(k+1+chnln1p)
             if(k.lt.chnln2) sz9y=old_rygenpep(k+3*chnln1p)+d2-old_rygenpep(k+1)
         			do d3=-move2(k),move2(k),.005d0
           				rcaz=old_rzgenpep(k+3*chnln1p)+d3-old_rzgenpep(k)
           				rnhz=old_rzgenpep(k+3*chnln1p)+d3-old_rzgenpep(k+chnln1p)
           				rcoz=old_rzgenpep(k+3*chnln1p)+d3-old_rzgenpep(k+2*chnln1p)
             if(k.ge.2) sz6z=old_rzgenpep(k+3*chnln1p)+d3-old_rzgenpep(k-1+2*chnln1p)
             if(k.ge.2) sz8z=old_rzgenpep(k+3*chnln1p)+d3-old_rzgenpep(k-1)
             if(k.lt.chnln2) sz7z=old_rzgenpep(k+3*chnln1p)+d3-old_rzgenpep(k+1+chnln1p)
             if(k.lt.chnln2) sz9z=old_rzgenpep(k+3*chnln1p)+d3-old_rzgenpep(k+1)
           				rca=dsqrt(rcax**2+rcay**2+rcaz**2)
           				rnh=dsqrt(rnhx**2+rnhy**2+rnhz**2)
           				rco=dsqrt(rcox**2+rcoy**2+rcoz**2)
             sz6=dsqrt(sz6x**2+sz6y**2+sz6z**2)
             sz7=dsqrt(sz7x**2+sz7y**2+sz7z**2)
             sz8=dsqrt(sz8x**2+sz8y**2+sz8z**2)
             sz9=dsqrt(sz9x**2+sz9y**2+sz9z**2)
!            sz10=dsqrt(sz10x**2+sz10y**2+sz10z**2)
                ds1 = dabs(bdln2(k)-rca)
                ds2 = dabs(bl_rn2(k)-rnh)
                ds3 = dabs(bl_rc2(k)-rco)
           				r3=dabs(bdln2(k)-rca)+dabs(bl_rn2(k)-rnh)+dabs(bl_rc2(k)-rco)
						 isz=0
             if(k.ge.2) then
             if(sz6.lt.sqz610cf(2,identity2(3*chnln2+k)-8)) isz=1
             if(sz8.lt.sqz610cf(1,identity2(3*chnln2+k)-8)) isz=1
             endif
             if(k.lt.chnln2) then
             if(sz7.lt.sqz610cf(3,identity2(3*chnln2+k)-8)) isz=1
             if(sz9.lt.sqz610cf(4,identity2(3*chnln2+k)-8)) isz=1
             endif
!            if(k.ge.3) then
!            if(sz10.lt.sqz610cf(5,identity1(3*chnln2+k)-8)) isz=1
!            endif
             if(isz.eq.0) id=id+1
             if((r3.lt.rmin).and.(isz.eq.0)) then
            					d1min=d1
          					    d2min=d2
               					d3min=d3
               					rmin=r3
           				endif
         			end do        
         		end do        
      end do        
		
       		old_rxgenpep(k+3*chnln1p)=old_rxgenpep(k+3*chnln1p)+d1min
       		old_rygenpep(k+3*chnln1p)=old_rygenpep(k+3*chnln1p)+d2min
       		old_rzgenpep(k+3*chnln1p)=old_rzgenpep(k+3*chnln1p)+d3min
       		rcax=old_rxgenpep(k+3*chnln1p)-old_rxgenpep(k)
       		rnhx=old_rxgenpep(k+3*chnln1p)-old_rxgenpep(k+chnln1p)
       		rcox=old_rxgenpep(k+3*chnln1p)-old_rxgenpep(k+2*chnln1p)
       		rcay=old_rygenpep(k+3*chnln1p)-old_rygenpep(k)
       		rnhy=old_rygenpep(k+3*chnln1p)-old_rygenpep(k+chnln1p)
       		rcoy=old_rygenpep(k+3*chnln1p)-old_rygenpep(k+2*chnln1p)
       		rcaz=old_rzgenpep(k+3*chnln1p)-old_rzgenpep(k)
       		rnhz=old_rzgenpep(k+3*chnln1p)-old_rzgenpep(k+chnln1p)
       		rcoz=old_rzgenpep(k+3*chnln1p)-old_rzgenpep(k+2*chnln1p)
       		rca=dsqrt(rcax**2+rcay**2+rcaz**2)
       		rnh=dsqrt(rnhx**2+rnhy**2+rnhz**2)
       		rco=dsqrt(rcox**2+rcoy**2+rcoz**2)
             if(k.ge.2) sz6x=old_rxgenpep(k+3*chnln2)-old_rxgenpep(k-1+2*chnln2)
             if(k.ge.2) sz8x=old_rxgenpep(k+3*chnln2)-old_rxgenpep(k-1)
             if(k.lt.chnln2) sz7x=old_rxgenpep(k+3*chnln2)-old_rxgenpep(k+1+chnln2)
             if(k.lt.chnln2) sz9x=old_rxgenpep(k+3*chnln2)-old_rxgenpep(k+1)
             if(k.ge.3) sz10x=old_rxgenpep(k+3*chnln2)-old_rxgenpep(k-2+2*chnln2)
             if(k.ge.2) sz6y=old_rygenpep(k+3*chnln2)-old_rygenpep(k-1+2*chnln2)
             if(k.ge.2) sz8y=old_rygenpep(k+3*chnln2)-old_rygenpep(k-1)
             if(k.lt.chnln2) sz7y=old_rygenpep(k+3*chnln2)-old_rygenpep(k+1+chnln2)
             if(k.lt.chnln2) sz9y=old_rygenpep(k+3*chnln2)-old_rygenpep(k+1)
             if(k.ge.3) sz10y=old_rygenpep(k+3*chnln2)-old_rygenpep(k-2+2*chnln2)
             if(k.ge.2) sz6z=old_rzgenpep(k+3*chnln2)-old_rzgenpep(k-1+2*chnln2)
             if(k.ge.2) sz8z=old_rzgenpep(k+3*chnln2)-old_rzgenpep(k-1)
             if(k.lt.chnln2) sz7z=old_rzgenpep(k+3*chnln2)-old_rzgenpep(k+1+chnln2)
             if(k.lt.chnln2) sz9z=old_rzgenpep(k+3*chnln2)-old_rzgenpep(k+1)
             if(k.ge.3) sz10z=old_rzgenpep(k+3*chnln2)-old_rzgenpep(k-2+2*chnln2)
             sz6=dsqrt(sz6x**2+sz6y**2+sz6z**2)
             sz7=dsqrt(sz7x**2+sz7y**2+sz7z**2)
             sz8=dsqrt(sz8x**2+sz8y**2+sz8z**2)
             sz9=dsqrt(sz9x**2+sz9y**2+sz9z**2)
             sz10=dsqrt(sz10x**2+sz10y**2+sz10z**2)
		write(6,'(i3,7(f8.4,1x))') k,rmin,d1min,d2min,d3min,rca,rnh,rco
		if (rmin .ge. 0.01) then
			write(6,*)'error in residue',k
			stop 'Configurating generation is terminated'
	        endif			
      end do

      pi = 4.d0*datan(1.d0)
      coll = 0
      t=0.d0
  
	!this gets rid of glycine residues
        	l=0
    	 do i=0,3
    		do k=1,chnln2
        		if(identity2(i*chnln2+k).ne.9) then
         			l=l+1
         			new_rx2(l)=old_rxgenpep(i*chnln1p+k)
         			new_ry2(l)=old_rygenpep(i*chnln1p+k)
         			new_rz2(l)=old_rzgenpep(i*chnln1p+k)
              endif
       		enddo
  	    enddo

      !open(7,file='inputs/chninfo-n2.data',status='unknown',form='formatted')
      do i=1,numbeads2
	      write(inpinfon2,*) new_rx2(i),new_ry2(i),new_rz2(i) ! xyz is stored in 4 groups
      enddo
      !close(unit=7)
	  
	!this makes a pdb to show the structure of one of each species
      !open(7,file='checks/configone2.pdb')
      write(checkcf2,*) numbeads2
      do i=1,numbeads2
       write(checkcf2,'(a6,i5,a3,1x,15x,3f8.3)') 'atom  ',i,'n', new_rx2(i), new_ry2(i), new_rz2(i)
      enddo
      !close(7)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!++++ READ single chain position information ++++++++++++++++++++!	  
999 if((exist_flag1).and.(exist_flag2)) write(6,*) 'reading chninfo.data and skip sidechain adjustment'
	rewind(inpinfon1)
      do i=1,numbeads1
	      read(inpinfon1,*) new_rxf(i),new_ryf(i),new_rzf(i)
      enddo	  
		
	 j = 0
          do i = 1,numbeads1
     			if ((i.ge.1) .and.(i.le.chnln1)) then
				j = j + 1
			       	new_rx(j)=new_rxf(i)
				new_ry(j)=new_ryf(i)
				new_rz(j)=new_rzf(i)
			elseif ((i.ge.(1+chnln1)) .and.(i.le.(chnln1+chnln1))) then 
				j = j + 1
			       	new_rx(j)=new_rxf(i)
				new_ry(j)=new_ryf(i)
				new_rz(j)=new_rzf(i)
   			elseif ((i.ge.(1+chnln1*2)) .and.(i.le.(chnln1+chnln1*2))) then 
				j = j + 1
			       	new_rx(j)=new_rxf(i)
				new_ry(j)=new_ryf(i)
				new_rz(j)=new_rzf(i)
   			elseif ((i.ge.(1+chnln1*3)) .and.(i.le.(chnln1+chnln1*3))) then 
				j = j + 1
			       	new_rx(j)=new_rxf(i)
				new_ry(j)=new_ryf(i)
				new_rz(j)=new_rzf(i)
               		endif
	enddo

	rewind(inpinfon2) 
      do i=1,numbeads2
	      read(inpinfon2,*) new_rxf2(i),new_ryf2(i),new_rzf2(i)
      enddo
	 j = 0
          do i = 1,numbeads2
     			if ((i.ge.1) .and.(i.le.chnln2)) then
				j = j + 1
			       	new_rx2(j)=new_rxf2(i)
				new_ry2(j)=new_ryf2(i)
				new_rz2(j)=new_rzf2(i)
			elseif ((i.ge.(1+chnln2)) .and.(i.le.(chnln2+chnln2))) then 
				j = j + 1
			        new_rx2(j)=new_rxf2(i)
				new_ry2(j)=new_ryf2(i)
				new_rz2(j)=new_rzf2(i)
   			elseif ((i.ge.(1+chnln2*2)) .and.(i.le.(chnln2+chnln2*2))) then 
				j = j + 1
			       	new_rx2(j)=new_rxf2(i)
				new_ry2(j)=new_ryf2(i)
				new_rz2(j)=new_rzf2(i)
   			elseif ((i.ge.(1+chnln2*3)) .and.(i.le.(chnln2+chnln2*3))) then 
				j = j + 1
			       	new_rx2(j)=new_rxf2(i)
				new_ry2(j)=new_ryf2(i)
				new_rz2(j)=new_rzf2(i)
               		endif
	enddo	  
	!this gets the cooridantes for the peptide relative to the first bead
      do i=1, numbeads1-1
      	do j=i+1, numbeads1
		dislocx(j) = new_rx(j) - new_rx(i)
		dislocy(j) = new_ry(j) - new_ry(i)
		dislocz(j) = new_rz(j) - new_rz(i)
      	enddo
      enddo
      do i=1,numbeads1
		xtmp(i)=new_rx(i)
		ytmp(i)=new_ry(i)
		ztmp(i)=new_rz(i)
      enddo
      do i=1, numbeads2-1
      	do j=i+1, numbeads2
		dislocx2(j) = new_rx2(j) - new_rx2(i)
		dislocy2(j) = new_ry2(j) - new_ry2(i)
		dislocz2(j) = new_rz2(j) - new_rz2(i)
      	enddo
      enddo
      do i=1,numbeads2
		xtmp(nop1+i)=new_rx2(i)
		ytmp(nop1+i)=new_ry2(i)
		ztmp(nop1+i)=new_rz2(i)
      enddo
      	!Read well width and well depth for all 20aa from the forcefield
      	do i=1,400
      		read(ffha55a,711) iiii,jjjj,bdtemp,wltemp
711  		format(2(2x,i2,2x),2(f6.3,2x))
       		bdsgencf(iiii,jjjj)=bdtemp
       		welcf(iiii,jjjj)=wltemp
      	end do

      if(bdsgencf(29,29)/2 .gt. 2.000) then
		sig_max=2.000+bdsgencf(29,29)/2
      else
		sig_max=4.0D0
      endif
      write(6,*)'Maximum Sigma',sig_max
      do i=1,29
        do j=1,29
			sig_ij=bdsgencf(i,j)
	       	if (sig_ij .gt. sig_max) then
                  		sig_max = sig_ij
      write(6,*)'Maximum Sigma',sig_max
              	endif
      enddo
      enddo
		
	!this selects a random position, build the peptide in that location, and checks for overlaps
	!if there is an overlap it repeats
      do i=1,nc
993		xtmp((i-1)*numbeads1+1)=(drandm(0))*boxl
		xtmp((i-1)*numbeads1+1)=xtmp((i-1)*numbeads1+1)-boxl*ANINT(xtmp((i-1)*numbeads1+1)/boxl)
		ytmp((i-1)*numbeads1+1)=(drandm(0))*boxl
		ytmp((i-1)*numbeads1+1)=ytmp((i-1)*numbeads1+1)-boxl*ANINT(ytmp((i-1)*numbeads1+1)/boxl)
		ztmp((i-1)*numbeads1+1)=(drandm(0))*boxl
		ztmp((i-1)*numbeads1+1)=ztmp((i-1)*numbeads1+1)-boxl*ANINT(ztmp((i-1)*numbeads1+1)/boxl)
       		do j=2,numbeads1
			xtmp((i-1)*numbeads1+j)=xtmp((i-1)*numbeads1+j-1)+dislocx(j)
			xtmp((i-1)*numbeads1+j)=xtmp((i-1)*numbeads1+j)-boxl*ANINT(xtmp((i-1)*numbeads1+j)/boxl)
			ytmp((i-1)*numbeads1+j)=ytmp((i-1)*numbeads1+j-1)+dislocy(j)
			ytmp((i-1)*numbeads1+j)=ytmp((i-1)*numbeads1+j)-boxl*ANINT(ytmp((i-1)*numbeads1+j)/boxl)
			ztmp((i-1)*numbeads1+j)=ztmp((i-1)*numbeads1+j-1)+dislocz(j)
			ztmp((i-1)*numbeads1+j)=ztmp((i-1)*numbeads1+j)-boxl*ANINT(ztmp((i-1)*numbeads1+j)/boxl)
       enddo
        do k=(i-1)*numbeads1+1,i*numbeads1
	      		do l=1,(i-1)*numbeads1
                       		rxij=xtmp(k)-xtmp(l)
                       		ryij=ytmp(k)-ytmp(l) 
                       		rzij=ztmp(k)-ztmp(l)
                      	 	rxij=rxij-boxl*aint((rxij/boxl)+.5)
                       		ryij=ryij-boxl*aint((ryij/boxl)+.5)
                       		rzij=rzij-boxl*aint((rzij/boxl)+.5)
	               		rijsq=rxij*rxij+ryij*ryij+rzij*rzij
	               		rijsq=rijsq*1.0000000001d0
                    if (rijsq .le. (5)**2) goto 993
           enddo
       enddo
      enddo
      do i=1,nc2
9993		xtmp(nop1+(i-1)*numbeads2+1)=(drandm(0))*boxl
		xtmp(nop1+(i-1)*numbeads2+1)=xtmp(nop1+(i-1)*numbeads2+1)-boxl*ANINT(xtmp(nop1+(i-1)*numbeads2+1)/boxl)
		ytmp(nop1+(i-1)*numbeads2+1)=(drandm(0))*boxl
		ytmp(nop1+(i-1)*numbeads2+1)=ytmp(nop1+(i-1)*numbeads2+1)-boxl*ANINT(ytmp(nop1+(i-1)*numbeads2+1)/boxl)
		ztmp(nop1+(i-1)*numbeads2+1)=(drandm(0))*boxl
		ztmp(nop1+(i-1)*numbeads2+1)=ztmp(nop1+(i-1)*numbeads2+1)-boxl*ANINT(ztmp(nop1+(i-1)*numbeads2+1)/boxl)
       		do j=2,numbeads2
			xtmp(nop1+(i-1)*numbeads2+j)=xtmp(nop1+(i-1)*numbeads2+j-1)+dislocx2(j)
			xtmp(nop1+(i-1)*numbeads2+j)=xtmp(nop1+(i-1)*numbeads2+j)-boxl*ANINT(xtmp(nop1+(i-1)*numbeads2+j)/boxl)
			ytmp(nop1+(i-1)*numbeads2+j)=ytmp(nop1+(i-1)*numbeads2+j-1)+dislocy2(j)
			ytmp(nop1+(i-1)*numbeads2+j)=ytmp(nop1+(i-1)*numbeads2+j)-boxl*ANINT(ytmp(nop1+(i-1)*numbeads2+j)/boxl)
			ztmp(nop1+(i-1)*numbeads2+j)=ztmp(nop1+(i-1)*numbeads2+j-1)+dislocz2(j)
			ztmp(nop1+(i-1)*numbeads2+j)=ztmp(nop1+(i-1)*numbeads2+j)-boxl*ANINT(ztmp(nop1+(i-1)*numbeads2+j)/boxl)
       enddo
        do k=(i-1)*numbeads2+1,i*numbeads2
	      		do l=1,(i-1)*numbeads2
                       		rxij=xtmp(nop1+k)-xtmp(nop1+l)
                       		ryij=ytmp(nop1+k)-ytmp(nop1+l) 
                       		rzij=ztmp(nop1+k)-ztmp(nop1+l)
                      	 	rxij=rxij-boxl*aint((rxij/boxl)+.5)
                       		ryij=ryij-boxl*aint((ryij/boxl)+.5)
                       		rzij=rzij-boxl*aint((rzij/boxl)+.5)
	               		rijsq=rxij*rxij+ryij*ryij+rzij*rzij
	               		rijsq=rijsq*1.0000000001d0
                    if (rijsq .le. (5)**2) goto 9993
           enddo
       enddo
      enddo
	!write full pdb file to verify configuration
      do i=1,nop1+nop2
       write(checkcf,'(a6,i5,a3,1x,15x,3f8.3)') 'ATOM  ',i,'N', xtmp(i), ytmp(i), ztmp(i)
      enddo
      !close(7)

!      open(7,file='checks/emblem.input')
!	j=0
!      do i=1,chnln1
!	j=j+1
!		write(7,'(2i5,15x,3f8.3,i5)') j,j, xtmp(chnln1+i), ytmp(chnln1+i), ztmp(chnln1+i),identity1(chnln1+i)
!	j=j+1
!		write(7,'(2i5,15x,3f8.3,i5)') j,j, xtmp(i), ytmp(i), ztmp(i), identity1(i)
!	j=j+1
!		write(7,'(2i5,15x,3f8.3,i5)') j,j, xtmp(chnln1*2+i), ytmp(chnln1*2+i), ztmp(chnln1*2+i), identity1(chnln1*2+i)
!	j=j+1
!		write(7,'(2i5,15x,3f8.3,i5)') j,j, xtmp(chnln1*3+i), ytmp(chnln1*3+i), ztmp(chnln1*3+i), identity1(chnln1*3+i)
!      enddo
!         close(7)
   	
      do i=1,nop1+nop2
   		if(xtmp(i).lt.xmin) xmin=xtmp(i)
   		if(xtmp(i).gt.xmax) xmax=xtmp(i)
   		if(ytmp(i).lt.ymin) ymin=ytmp(i)
   		if(ytmp(i).gt.ymax) ymax=ytmp(i)
   		if(ztmp(i).lt.zmin) zmin=ztmp(i)
   		if(ztmp(i).gt.zmax) zmax=ztmp(i)
      end do
	
	!check is total bounds are within box length	
      write(6,*) 'These should be less than',boxl
      write(6,*) 'X-direction',abs(xmin)+abs(xmax)
      write(6,*) 'Y-direction',abs(ymin)+abs(ymax)
      write(6,*) 'Z-direction',abs(zmin)+abs(zmax)

	!write a config check
      write(checkcfout,*) coll,t
      write(checkcfout,*)xtmp/boxl,ytmp/boxl,ztmp/boxl
      !close(7)	
	
	!this writes the actual config input
        write(runcf) coll,t,xtmp/boxl,ytmp/boxl,ztmp/boxl
      !close(7)
	
	!read in masses and adjust for bead identity
	
	bmass(3)=bmass(20)
	
      do i=1,4
		bmass(i+4)=bmass(i)
      end do

      !close(7)

      do i=1,nopwg1
		bm1(i)=bmass(identity1(i))
      end do
      do i=1,nopwg2
		 bm2(i)=bmass(identity2(i))
      end do	  

	!write to make sure masses are correct
	!call chdir(rundir)
      !open(7,file='checks/masses.out')
       do i=1,nopwg1
        write(checkmass,*) i,identity1(i),bm1(i)
      enddo
       do i=1,nopwg2
        write(checkmass,*) nopwg1+i,identity2(i),bm2(i)
      enddo
       !close(7)

      write(6,*)'Desired reduced temperature',setemp/12

	!assign random velocities
      do k=1,nopwg1
       if(identity1(k).ne.9) then
         	x1=drandm(0)*0.75d0+0.2d0
         	x2=drandm(0)*0.75d0+0.2d0
         	y1=drandm(0)*0.75d0+0.2d0
         	y2=drandm(0)*0.75d0+0.2d0
         	svgcf(4,k)=dsqrt(-2.d0*log(x1))*cos(2.d0*pi*y1)
         	svgcf(5,k)=dsqrt(-2.d0*log(x1))*sin(2.d0*pi*y1)
         	svgcf(6,k)=dsqrt(-2.d0*log(x2))*cos(2.d0*pi*y2)
       endif
      end do
      do k=1,nopwg2
       if(identity2(k).ne.9) then
         	x1=drandm(0)*0.75d0+0.2d0
         	x2=drandm(0)*0.75d0+0.2d0
         	y1=drandm(0)*0.75d0+0.2d0
         	y2=drandm(0)*0.75d0+0.2d0
         	svgcf(4,nopwg1+k)=dsqrt(-2.d0*log(x1))*cos(2.d0*pi*y1)
         	svgcf(5,nopwg1+k)=dsqrt(-2.d0*log(x1))*sin(2.d0*pi*y1)
         	svgcf(6,nopwg1++k)=dsqrt(-2.d0*log(x2))*cos(2.d0*pi*y2)
       endif
      end do

	!scale velocities so that linear momentum = 0
      do k=1,nopwg1
       if(identity1(k).ne.9) then
         	sumx=sumx+svgcf(4,k)
         	sumy=sumy+svgcf(5,k)
         	sumz=sumz+svgcf(6,k)
      		endif
      end do
      do k=1,nopwg2
       if(identity2(k).ne.9) then
         	sumx=sumx+svgcf(4,nopwg1+k)
         	sumy=sumy+svgcf(5,nopwg1+k)
         	sumz=sumz+svgcf(6,nopwg1+k)
      		endif
      end do	
	
      do k=1,nopwg1
       if(identity1(k).ne.9) then
         	svgcf(4,k)=svgcf(4,k)-sumx/real(nop1+nop2)
         	svgcf(5,k)=svgcf(5,k)-sumy/real(nop1+nop2)
         	svgcf(6,k)=svgcf(6,k)-sumz/real(nop1+nop2)
     	endif
      end do
      do k=1,nopwg2
       if(identity2(k).ne.9) then
         	svgcf(4,nopwg1+k)=svgcf(4,nopwg1+k)-sumx/real(nop1+nop2)
         	svgcf(5,nopwg1+k)=svgcf(5,nopwg1+k)-sumy/real(nop1+nop2)
         	svgcf(6,nopwg1+k)=svgcf(6,nopwg1+k)-sumz/real(nop1+nop2)
     	endif
      end do	
	sumvel=0.d0
	
	!check summed velocities
      !open(7,file='checks/sumvelcheck.out')
      do k=1,nopwg1
       if(identity1(k).ne.9) then
			sumvel=sumvel+(svgcf(4,k)**2+svgcf(5,k)**2+svgcf(6,k)**2)/bm1(k)
       endif
      enddo  
      do k=1,nopwg2
       if(identity2(k).ne.9) then
			sumvel=sumvel+(svgcf(4,nopwg1+k)**2+svgcf(5,nopwg1+k)**2+svgcf(6,nopwg1+k)**2)/bm2(k)
       endif
      enddo 
	!here svgcf is actually momentum not velocity
	 
	tred=sumvel/3.d0/dble(nop1+nop2)
	const=dsqrt(setemp/tred)

       do k=1,nopwg1
       if(identity1(k).ne.9) then
     		svgcf(4,k)=const*svgcf(4,k)/bm1(k)
     		svgcf(5,k)=const*svgcf(5,k)/bm1(k)
     		svgcf(6,k)=const*svgcf(6,k)/bm1(k)
        endif
      enddo
       do k=1,nopwg2
       if(identity2(k).ne.9) then
     		svgcf(4,nopwg1+k)=const*svgcf(4,nopwg1+k)/bm2(k)
     		svgcf(5,nopwg1+k)=const*svgcf(5,nopwg1+k)/bm2(k)
     		svgcf(6,nopwg1+k)=const*svgcf(6,nopwg1+k)/bm2(k)
        endif
      enddo
  	
	sumvel=0.d0
      
      do k=1,nopwg1
        if(identity1(k).ne.9) then
			sumvel=sumvel+bm1(k)*(svgcf(4,k)**2+svgcf(5,k)**2+svgcf(6,k)**2)
			write(checksumvel,*)k,sumvel
      	endif
      enddo
      do k=1,nopwg2
        if(identity2(k).ne.9) then
			sumvel=sumvel+bm2(k)*(svgcf(4,nopwg1+k)**2+svgcf(5,nopwg1+k)**2+svgcf(6,nopwg1+k)**2)
			write(checksumvel,*)nopwg1+k,sumvel
      	endif
      enddo
      !close(7)
	
	tred=sumvel/3.d0/dble(nop1+nop2)
      write(6,*)'New reduced temperature',tred/12.0

	!this is the velocity input to the simulations
      !open(7,file='results/run0000.lastvel',status='unknown',form='unformatted') 
	l=0
      do k=1,nopwg1   
         if(identity1(k).ne.9) then
         	l=l+1
         	xtmp(l)=svgcf(4,k)
         	ytmp(l)=svgcf(5,k)
         	ztmp(l)=svgcf(6,k)
        endif
      enddo
      do k=1,nopwg2   
         if(identity2(k).ne.9) then
         	l=l+1
         	xtmp(l)=svgcf(4,nopwg1+k)
         	ytmp(l)=svgcf(5,nopwg1+k)
         	ztmp(l)=svgcf(6,nopwg1+k)
        endif
      enddo
      write(runlasvel) coll,xtmp,ytmp,ztmp

	!check velocities
	l=0
       do k=1,nopwg1  
       		if(identity1(k).ne.9) then
         		l=l+1
         		xtmp(l)=svgcf(4,k)
         		ytmp(l)=svgcf(5,k)
         		ztmp(l)=svgcf(6,k)
        		write(checkvel,*) k,identity1(k),xtmp(k),ytmp(k),ztmp(k)
    		endif
      enddo
       do k=1,nopwg2   
       		if(identity2(k).ne.9) then
         		l=l+1
         		xtmp(l)=svgcf(4,nopwg1+k)
         		ytmp(l)=svgcf(5,nopwg1+k)
         		ztmp(l)=svgcf(6,nopwg1+k)
        		write(checkvel,*) nopwg1+k,identity2(k),xtmp(nopwg1+k),ytmp(nopwg1+k),ztmp(nopwg1+k)
    		endif
      enddo
      

	!check HB parnters
      do i=1,numbeads1
		write(checkhb,*)i,bptnr1(i)
      enddo
      do i=1,numbeads2
		write(checkhb,*)nopwg1+i,bptnr2(i)
      enddo	  
	call gencf_filecloser()
      end
#include "files_opn_gencf.f"
#include "files_close_gencf.f"
