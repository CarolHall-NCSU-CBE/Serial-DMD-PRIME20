	use global
     	implicit none
      	integer, parameter :: k12=selected_int_kind(12)
        integer nb1p,nb21p,chnln1p
	parameter(nb1p=124)       ! number of beads 124, we don't change lines 6-8
	parameter(nb21p=124)       ! number of beads 124	  
        parameter(chnln1p=31)     ! chain length  31  
      integer(kind=k12) coll
      real*8 drandm,xr
      integer iflag
      external drandm,srand

      integer aa(nb),identity(nop),bptnr(nop2)
      integer aa2(nb3),identity2(nop3),bptnr2(nop4)	  
      integer i,j,k,l,iii,iiii,numbeads,jjjj

      real*8 old_rx(nb1p), old_ry(nb1p), old_rz(nb1p)     
      real*8 dave_rx(124), dave_ry(124), dave_rz(124)     	  
      real*8 new_rx(nb2), new_ry(nb2), new_rz(nb2)
      real*8 new_rx2(nb4), new_ry2(nb4), new_rz2(nb4)	  
      real*8 new_rxf(nb2), new_ryf(nb2), new_rzf(nb2) 	  
      real*8 new_rxf2(nb4), new_ryf2(nb4), new_rzf2(nb4)
      real*8 xtmp(nop2+nop4),ytmp(nop2+nop4),ztmp(nop2+nop4)            
      real*8 sv(6,nop+nop3),boxl_orig,t
      real*8 xmin,xmax,ymin,ymax,zmin,zmax
      real*8 dislocx(nb2),dislocy(nb2),dislocz(nb2)
      real*8 dislocx2(nb4),dislocy2(nb4),dislocz2(nb4)
      real sig_ij,sig_max,bds(400,400),wel(400,400),bdtemp,wltemp
 
      real*8 rmin,r3,rca,rnh,rco,move(chnln),move2(chnln2),DS1, DS2, DS3                              
      real*8 rcax,rcay,rcaz,rnhx,rnhy,rnhz,rcox,rcoy,rcoz
      real*8 del,pi,setemp,reducedtemp
      real*8 drca(20),drnh(20),drco(20)
      real*8 del_rca(20),del_rnh(20),del_rco(20)
      real*8 bdln(chnln),bl_rn(chnln),bl_rc(chnln),del_bdln(chnln),del_blrn(chnln),del_blrc(chnln)
      real*8 bdln2(chnln),bl_rn2(chnln),bl_rc2(chnln),del_bdln2(chnln),del_blrn2(chnln),del_blrc2(chnln)
      real*8 scalerca,rx2,ry2,rz2
      real*8 d1,d2,d3,d1min,d2min,d3min
      real*8 rxij,ryij,rzij,rijsq,boxl
      REAL*8 SQZ610(5,20)
      REAL*8 SZ6,SZ7,SZ8,SZ9,SZ10
      REAL*8 SZ6X,SZ7X,SZ8X,SZ9X,SZ10X
      REAL*8 SZ6Y,SZ7Y,SZ8Y,SZ9Y,SZ10Y
      REAL*8 SZ6Z,SZ7Z,SZ8Z,SZ9Z,SZ10Z
      INTEGER ISZ,ID
	  
      real*8 bm(nop),bm2(nop3),bmass_temp,bmass(28)                      
      real*8 sumvel,tred,const,sumx,sumy,sumz,x1,x2,y1,y2

!	To creat input file
	char*31 seq1, seq2       

      logical exist_flag1,exist_flag2
      inquire(file = 'chninfo-n1.data', exist = exist_flag1)	 
      inquire(file = 'chninfo-n2.data', exist = exist_flag2)	  
	  
      iflag = 1331171207 !1058472402   !1970711301   !1331171207
      xr=drandm(iflag)
      call srand(iflag)
      boxl=159.0D0
    
      reducedtemp=0.50
      setemp = reducedtemp*12  
	  
      open(7,file='inputs/peptidex.inp',status='old')  
      read(7,*) numbeads,boxl_orig,old_rx
      close(7)
      open(7,file='inputs/peptidey.inp',status='old')  
      read(7,*) numbeads,boxl_orig,old_ry
      close(7)
      open(7,file='inputs/peptidez.inp',status='old')  
      read(7,*) numbeads,boxl_orig,old_rz
      close(7)
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
!      do i=1,chnln
!      read(7,'(30x,3(f8.3))') old_rx(i+chnln),old_ry(i+chnln),old_rz(i+chnln)
!      read(7,'(30x,3(f8.3))') old_rx(i),old_ry(i),old_rz(i)
!      read(7,'(30x,3(f8.3))') old_rx(i+2*chnln),old_ry(i+2*chnln),old_rz(i+2*chnln)
!      read(7,'(30x,3(f8.3))') old_rx(i+3*chnln),old_ry(i+3*chnln),old_rz(i+3*chnln)
!      end do
!      close(7)

	!read in identities
      open(7,file='inputs/identity.inp',status='old')   
      read (7,*) aa
      close(7)
      
      do l=1,nop,nb
       do k=1,nb 
		identity(l+k-1)=aa(k)
      enddo
      enddo

      open(7,file='inputs/identity2.inp',status='old')   
      read (7,*) aa2
      close(7)
      
      do l=1,nop3,nb3
       do k=1,nb3 
		identity2(l+k-1)=aa2(k)
      enddo
      enddo	  
	  
	!write file to verify identities
      open(7,file='checks/identity.out')
      do i=1,nop
		write(7,*) i,identity(i)
      enddo
      do i=1,nop3
		write(7,*) i,identity2(i)
      enddo	  
      close(7)
	
	!this writes the hydrophobicity inputs
      open(7,file='parameters/hp1.inp')
      do i=1,nb
      if (identity(i) .le. 4) then
			write(7,*)'0'
      elseif (identity(i) .gt. 9) then
			write(7,*)'1'			
      endif
      enddo
      close(7)
      open(7,file='parameters/hp2.inp')
      do i=1,nb3
      if (identity2(i) .le. 4) then
			write(7,*)'0'
      elseif (identity2(i) .gt. 9) then
			write(7,*)'1'			
      endif
      enddo
      close(7)

	!this writes the sidechain inputs
      open(7,file='parameters/firstside1.data')
      do i=chnln*3+1,nb
      if (identity(i) .eq. 9) then
			write(7,*)'0'
      else
			write(7,*)'1'
      endif
      enddo
      close(7)
      open(7,file='parameters/firstside2.data')
      do i=chnln2*3+1,nb3
      if (identity2(i) .eq. 9) then
			write(7,*)'0'
      else
			write(7,*)'1'
      endif
      enddo
      close(7)
	  
	!this writes the identity input for the simulation
      open(7,file='parameters/identity.inp')   
      do l=1,nb
      	if (identity(l) .ne. 9) write(7,*)identity(l)
      enddo
      do l=1,nb3
      	if (identity2(l) .ne. 9) write(7,*)identity2(l)
      enddo
      close(7)

	del = 0.02375d0

	!read in l/d-isomer pseudobond constraints
      open(7,file='parameters/rcarnrco.data',status='old')
      do i=1,20
		read(7,"(6(f6.3,2x))") drca(i),drnh(i),drco(i),del_rca(i),del_rnh(i),del_rco(i)
	      if(del_rca(i).lt.del) del_rca(i)=del
      		if(del_rnh(i).lt.del) del_rnh(i)=del
      		if(del_rco(i).lt.del) del_rco(i)=del
      end do
      close(7)

	!reassign pseudobond identities
      do i=1,chnln
		iii=identity(3*chnln+i)-8
        	bdln(i)=drca(iii)
        	bl_rn(i)=drnh(iii)
        	bl_rc(i)=drco(iii)
        	del_bdln(i)=del_rca(iii)
        	del_blrn(i)=del_rnh(iii)
        	del_blrc(i)=del_rco(iii)
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

       OPEN(UNIT=7,FILE='parameters/sqz6to10.data',STATUS='OLD')
      DO I=1,20
       READ(7,'(5(f6.3,2x))') (SQZ610(K,I),K=1,5)
      END DO
       CLOSE(UNIT=7)
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if((exist_flag1).and.(exist_flag2)) goto 999	! Yiming
	
        do i=1,chnln
           move(i) = 3.5 
            enddo
	
         !this entire section takes each sidechain and moves the location of the bead in a box until it satisfies the pseudobond contraints
      do k=3*chnln1p+1,3*chnln1p+chnln
		rx2=(old_rx(k)-old_rx(k-3*chnln1p))
		ry2=(old_ry(k)-old_ry(k-3*chnln1p))
		rz2=(old_rz(k)-old_rz(k-3*chnln1p))
		scalerca=bdln(k-3*chnln1p)/dsqrt(rx2**2+ry2**2+rz2**2)
      		old_rx(k)=old_rx(k-3*chnln1p)+scalerca*rx2
      		old_ry(k)=old_ry(k-3*chnln1p)+scalerca*ry2
      		old_rz(k)=old_rz(k-3*chnln1p)+scalerca*rz2
      end do
                ds1 = 0.02
                ds2 = 0.02
                ds3 = 0.02
      write(6,*)
      write(6,'(a50)') "first value should be less than the 0.01 tolerance"
      write(6,'(a62)') "if not, adjust d1 d2 d3 values in sidechain adjustment section"

      do k=1,chnln
		rmin=100.d0		
         	do d1=-move(k),move(k),0.005d0
			rcax=old_rx(k+3*chnln1p)+d1-old_rx(k)
          		rnhx=old_rx(k+3*chnln1p)+d1-old_rx(k+chnln1p)
          		rcox=old_rx(k+3*chnln1p)+d1-old_rx(k+2*chnln1p)
             if(k.ge.2) sz6x=old_rx(k+3*chnln1p)+d1-old_rx(k-1+2*chnln1p)
             if(k.ge.2) sz8x=old_rx(k+3*chnln1p)+d1-old_rx(k-1)
             if(k.lt.chnln) sz7x=old_rx(k+3*chnln1p)+d1-old_rx(k+1+chnln1p)
             if(k.lt.chnln) sz9x=old_rx(k+3*chnln1p)+d1-old_rx(k+1)
         		do d2=-move(k),move(k),.005d0
           			rcay=old_ry(k+3*chnln1p)+d2-old_ry(k)
           			rnhy=old_ry(k+3*chnln1p)+d2-old_ry(k+chnln1p)
           			rcoy=old_ry(k+3*chnln1p)+d2-old_ry(k+2*chnln1p)
             if(k.ge.2) sz6y=old_ry(k+3*chnln1p)+d2-old_ry(k-1+2*chnln1p)
             if(k.ge.2) sz8y=old_ry(k+3*chnln1p)+d2-old_ry(k-1)
             if(k.lt.chnln) sz7y=old_ry(k+3*chnln1p)+d2-old_ry(k+1+chnln1p)
             if(k.lt.chnln) sz9y=old_ry(k+3*chnln1p)+d2-old_ry(k+1)
         			do d3=-move(k),move(k),.005d0
           				rcaz=old_rz(k+3*chnln1p)+d3-old_rz(k)
           				rnhz=old_rz(k+3*chnln1p)+d3-old_rz(k+chnln1p)
           				rcoz=old_rz(k+3*chnln1p)+d3-old_rz(k+2*chnln1p)
             if(k.ge.2) sz6z=old_rz(k+3*chnln1p)+d3-old_rz(k-1+2*chnln1p)
             if(k.ge.2) sz8z=old_rz(k+3*chnln1p)+d3-old_rz(k-1)
             if(k.lt.chnln) sz7z=old_rz(k+3*chnln1p)+d3-old_rz(k+1+chnln1p)
             if(k.lt.chnln) sz9z=old_rz(k+3*chnln1p)+d3-old_rz(k+1)
           				rca=dsqrt(rcax**2+rcay**2+rcaz**2)
           				rnh=dsqrt(rnhx**2+rnhy**2+rnhz**2)
           				rco=dsqrt(rcox**2+rcoy**2+rcoz**2)
             sz6=dsqrt(sz6x**2+sz6y**2+sz6z**2)
             sz7=dsqrt(sz7x**2+sz7y**2+sz7z**2)
             sz8=dsqrt(sz8x**2+sz8y**2+sz8z**2)
             sz9=dsqrt(sz9x**2+sz9y**2+sz9z**2)
!            sz10=dsqrt(sz10x**2+sz10y**2+sz10z**2)
                ds1 = dabs(bdln(k)-rca)
                ds2 = dabs(bl_rn(k)-rnh)
                ds3 = dabs(bl_rc(k)-rco)
           				r3=dabs(bdln(k)-rca)+dabs(bl_rn(k)-rnh)+dabs(bl_rc(k)-rco)
						 isz=0
             if(k.ge.2) then
             if(sz6.lt.sqz610(2,identity(3*chnln+k)-8)) isz=1
             if(sz8.lt.sqz610(1,identity(3*chnln+k)-8)) isz=1
             endif
             if(k.lt.chnln) then
             if(sz7.lt.sqz610(3,identity(3*chnln+k)-8)) isz=1
             if(sz9.lt.sqz610(4,identity(3*chnln+k)-8)) isz=1
             endif
!            if(k.ge.3) then
!            if(sz10.lt.sqz610(5,identity(3*chnln+k)-8)) isz=1
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
		
       		old_rx(k+3*chnln1p)=old_rx(k+3*chnln1p)+d1min
       		old_ry(k+3*chnln1p)=old_ry(k+3*chnln1p)+d2min
       		old_rz(k+3*chnln1p)=old_rz(k+3*chnln1p)+d3min
       		rcax=old_rx(k+3*chnln1p)-old_rx(k)
       		rnhx=old_rx(k+3*chnln1p)-old_rx(k+chnln1p)
       		rcox=old_rx(k+3*chnln1p)-old_rx(k+2*chnln1p)
       		rcay=old_ry(k+3*chnln1p)-old_ry(k)
       		rnhy=old_ry(k+3*chnln1p)-old_ry(k+chnln1p)
       		rcoy=old_ry(k+3*chnln1p)-old_ry(k+2*chnln1p)
       		rcaz=old_rz(k+3*chnln1p)-old_rz(k)
       		rnhz=old_rz(k+3*chnln1p)-old_rz(k+chnln1p)
       		rcoz=old_rz(k+3*chnln1p)-old_rz(k+2*chnln1p)
       		rca=dsqrt(rcax**2+rcay**2+rcaz**2)
       		rnh=dsqrt(rnhx**2+rnhy**2+rnhz**2)
       		rco=dsqrt(rcox**2+rcoy**2+rcoz**2)
             if(k.ge.2) sz6x=old_rx(k+3*chnln)-old_rx(k-1+2*chnln)
             if(k.ge.2) sz8x=old_rx(k+3*chnln)-old_rx(k-1)
             if(k.lt.chnln) sz7x=old_rx(k+3*chnln)-old_rx(k+1+chnln)
             if(k.lt.chnln) sz9x=old_rx(k+3*chnln)-old_rx(k+1)
             if(k.ge.3) sz10x=old_rx(k+3*chnln)-old_rx(k-2+2*chnln)
             if(k.ge.2) sz6y=old_ry(k+3*chnln)-old_ry(k-1+2*chnln)
             if(k.ge.2) sz8y=old_ry(k+3*chnln)-old_ry(k-1)
             if(k.lt.chnln) sz7y=old_ry(k+3*chnln)-old_ry(k+1+chnln)
             if(k.lt.chnln) sz9y=old_ry(k+3*chnln)-old_ry(k+1)
             if(k.ge.3) sz10y=old_ry(k+3*chnln)-old_ry(k-2+2*chnln)
             if(k.ge.2) sz6z=old_rz(k+3*chnln)-old_rz(k-1+2*chnln)
             if(k.ge.2) sz8z=old_rz(k+3*chnln)-old_rz(k-1)
             if(k.lt.chnln) sz7z=old_rz(k+3*chnln)-old_rz(k+1+chnln)
             if(k.lt.chnln) sz9z=old_rz(k+3*chnln)-old_rz(k+1)
             if(k.ge.3) sz10z=old_rz(k+3*chnln)-old_rz(k-2+2*chnln)
             sz6=dsqrt(sz6x**2+sz6y**2+sz6z**2)
             sz7=dsqrt(sz7x**2+sz7y**2+sz7z**2)
             sz8=dsqrt(sz8x**2+sz8y**2+sz8z**2)
             sz9=dsqrt(sz9x**2+sz9y**2+sz9z**2)
             sz10=dsqrt(sz10x**2+sz10y**2+sz10z**2)			
            if (rmin .ge. 0.01) write(6,*)'error in residue',k
        	write(6,'(i3,7(f8.4,1x))') k,rmin,d1min,d2min,d3min,rca,rnh,rco
      end do

      pi = 4.d0*datan(1.d0)
      coll = 0
      t=0.d0
  
	!this gets rid of glycine residues
        	l=0
    	 do i=0,3
    		do k=1,chnln
        		if(identity(i*chnln+k).ne.9) then
         			l=l+1
         			new_rx(l)=old_rx(i*chnln1p+k)
         			new_ry(l)=old_ry(i*chnln1p+k)
         			new_rz(l)=old_rz(i*chnln1p+k)
              endif
       		enddo
  	    enddo
		
      open(7,file='chninfo-n1.data',status='unknown',form='formatted')
      do i=1,nb2
	      write(7,*) new_rx(i),new_ry(i),new_rz(i) ! xyz is stored in 4 groups
      enddo
      close(unit=7)
	  
	!this makes a pdb to show the structure of one of each species
      open(7,file='checks/configone1.pdb')
      write(7,*) nb2
      do i=1,nb2
       write(7,'(a6,i5,a3,1x,15x,3f8.3)') 'atom  ',i,'n', new_rx(i), new_ry(i), new_rz(i)
      enddo
      close(7)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=1,chnln2
           move2(i) = 2.5 
            enddo
	
    	d1min=0.0
        d2min=0.0
     	d3min=0.0
         rmin=0.0
	
      open(7,file='inputs/peptidex.inp',status='old')  
      read(7,*) numbeads,boxl_orig,old_rx
      close(7)
      open(7,file='inputs/peptidey.inp',status='old')  
      read(7,*) numbeads,boxl_orig,old_ry
      close(7)
      open(7,file='inputs/peptidez.inp',status='old')  
      read(7,*) numbeads,boxl_orig,old_rz
      close(7)
	
         !this entire section takes each sidechain and moves the location of the bead in a box until it satisfies the pseudobond contraints
      do k=3*chnln1p+1,3*chnln1p+chnln2
		rx2=(old_rx(k)-old_rx(k-3*chnln1p))
		ry2=(old_ry(k)-old_ry(k-3*chnln1p))
		rz2=(old_rz(k)-old_rz(k-3*chnln1p))
		scalerca=bdln2(k-3*chnln1p)/dsqrt(rx2**2+ry2**2+rz2**2)
      		old_rx(k)=old_rx(k-3*chnln1p)+scalerca*rx2
      		old_ry(k)=old_ry(k-3*chnln1p)+scalerca*ry2
      		old_rz(k)=old_rz(k-3*chnln1p)+scalerca*rz2
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
			rcax=old_rx(k+3*chnln1p)+d1-old_rx(k)
          		rnhx=old_rx(k+3*chnln1p)+d1-old_rx(k+chnln1p)
          		rcox=old_rx(k+3*chnln1p)+d1-old_rx(k+2*chnln1p)
             if(k.ge.2) sz6x=old_rx(k+3*chnln1p)+d1-old_rx(k-1+2*chnln1p)
             if(k.ge.2) sz8x=old_rx(k+3*chnln1p)+d1-old_rx(k-1)
             if(k.lt.chnln2) sz7x=old_rx(k+3*chnln1p)+d1-old_rx(k+1+chnln1p)
             if(k.lt.chnln2) sz9x=old_rx(k+3*chnln1p)+d1-old_rx(k+1)
         		do d2=-move2(k),move2(k),.005d0
           			rcay=old_ry(k+3*chnln1p)+d2-old_ry(k)
           			rnhy=old_ry(k+3*chnln1p)+d2-old_ry(k+chnln1p)
           			rcoy=old_ry(k+3*chnln1p)+d2-old_ry(k+2*chnln1p)
             if(k.ge.2) sz6y=old_ry(k+3*chnln1p)+d2-old_ry(k-1+2*chnln1p)
             if(k.ge.2) sz8y=old_ry(k+3*chnln1p)+d2-old_ry(k-1)
             if(k.lt.chnln2) sz7y=old_ry(k+3*chnln1p)+d2-old_ry(k+1+chnln1p)
             if(k.lt.chnln2) sz9y=old_ry(k+3*chnln1p)+d2-old_ry(k+1)
         			do d3=-move2(k),move2(k),.005d0
           				rcaz=old_rz(k+3*chnln1p)+d3-old_rz(k)
           				rnhz=old_rz(k+3*chnln1p)+d3-old_rz(k+chnln1p)
           				rcoz=old_rz(k+3*chnln1p)+d3-old_rz(k+2*chnln1p)
             if(k.ge.2) sz6z=old_rz(k+3*chnln1p)+d3-old_rz(k-1+2*chnln1p)
             if(k.ge.2) sz8z=old_rz(k+3*chnln1p)+d3-old_rz(k-1)
             if(k.lt.chnln2) sz7z=old_rz(k+3*chnln1p)+d3-old_rz(k+1+chnln1p)
             if(k.lt.chnln2) sz9z=old_rz(k+3*chnln1p)+d3-old_rz(k+1)
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
             if(sz6.lt.sqz610(2,identity2(3*chnln2+k)-8)) isz=1
             if(sz8.lt.sqz610(1,identity2(3*chnln2+k)-8)) isz=1
             endif
             if(k.lt.chnln2) then
             if(sz7.lt.sqz610(3,identity2(3*chnln2+k)-8)) isz=1
             if(sz9.lt.sqz610(4,identity2(3*chnln2+k)-8)) isz=1
             endif
!            if(k.ge.3) then
!            if(sz10.lt.sqz610(5,identity(3*chnln2+k)-8)) isz=1
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
		
       		old_rx(k+3*chnln1p)=old_rx(k+3*chnln1p)+d1min
       		old_ry(k+3*chnln1p)=old_ry(k+3*chnln1p)+d2min
       		old_rz(k+3*chnln1p)=old_rz(k+3*chnln1p)+d3min
       		rcax=old_rx(k+3*chnln1p)-old_rx(k)
       		rnhx=old_rx(k+3*chnln1p)-old_rx(k+chnln1p)
       		rcox=old_rx(k+3*chnln1p)-old_rx(k+2*chnln1p)
       		rcay=old_ry(k+3*chnln1p)-old_ry(k)
       		rnhy=old_ry(k+3*chnln1p)-old_ry(k+chnln1p)
       		rcoy=old_ry(k+3*chnln1p)-old_ry(k+2*chnln1p)
       		rcaz=old_rz(k+3*chnln1p)-old_rz(k)
       		rnhz=old_rz(k+3*chnln1p)-old_rz(k+chnln1p)
       		rcoz=old_rz(k+3*chnln1p)-old_rz(k+2*chnln1p)
       		rca=dsqrt(rcax**2+rcay**2+rcaz**2)
       		rnh=dsqrt(rnhx**2+rnhy**2+rnhz**2)
       		rco=dsqrt(rcox**2+rcoy**2+rcoz**2)
             if(k.ge.2) sz6x=old_rx(k+3*chnln2)-old_rx(k-1+2*chnln2)
             if(k.ge.2) sz8x=old_rx(k+3*chnln2)-old_rx(k-1)
             if(k.lt.chnln2) sz7x=old_rx(k+3*chnln2)-old_rx(k+1+chnln2)
             if(k.lt.chnln2) sz9x=old_rx(k+3*chnln2)-old_rx(k+1)
             if(k.ge.3) sz10x=old_rx(k+3*chnln2)-old_rx(k-2+2*chnln2)
             if(k.ge.2) sz6y=old_ry(k+3*chnln2)-old_ry(k-1+2*chnln2)
             if(k.ge.2) sz8y=old_ry(k+3*chnln2)-old_ry(k-1)
             if(k.lt.chnln2) sz7y=old_ry(k+3*chnln2)-old_ry(k+1+chnln2)
             if(k.lt.chnln2) sz9y=old_ry(k+3*chnln2)-old_ry(k+1)
             if(k.ge.3) sz10y=old_ry(k+3*chnln2)-old_ry(k-2+2*chnln2)
             if(k.ge.2) sz6z=old_rz(k+3*chnln2)-old_rz(k-1+2*chnln2)
             if(k.ge.2) sz8z=old_rz(k+3*chnln2)-old_rz(k-1)
             if(k.lt.chnln2) sz7z=old_rz(k+3*chnln2)-old_rz(k+1+chnln2)
             if(k.lt.chnln2) sz9z=old_rz(k+3*chnln2)-old_rz(k+1)
             if(k.ge.3) sz10z=old_rz(k+3*chnln2)-old_rz(k-2+2*chnln2)
             sz6=dsqrt(sz6x**2+sz6y**2+sz6z**2)
             sz7=dsqrt(sz7x**2+sz7y**2+sz7z**2)
             sz8=dsqrt(sz8x**2+sz8y**2+sz8z**2)
             sz9=dsqrt(sz9x**2+sz9y**2+sz9z**2)
             sz10=dsqrt(sz10x**2+sz10y**2+sz10z**2)			
            if (rmin .ge. 0.01) write(6,*)'error in residue',k
        	write(6,'(i3,7(f8.4,1x))') k,rmin,d1min,d2min,d3min,rca,rnh,rco
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
         			new_rx2(l)=old_rx(i*chnln1p+k)
         			new_ry2(l)=old_ry(i*chnln1p+k)
         			new_rz2(l)=old_rz(i*chnln1p+k)
              endif
       		enddo
  	    enddo
		
      open(7,file='chninfo-n2.data',status='unknown',form='formatted')
      do i=1,nb4
	      write(7,*) new_rx2(i),new_ry2(i),new_rz2(i) ! xyz is stored in 4 groups
      enddo
      close(unit=7)
	  
	!this makes a pdb to show the structure of one of each species
      open(7,file='checks/configone2.pdb')
      write(7,*) nb4
      do i=1,nb4
       write(7,'(a6,i5,a3,1x,15x,3f8.3)') 'atom  ',i,'n', new_rx2(i), new_ry2(i), new_rz2(i)
      enddo
      close(7)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!++++ READ single chain position information ++++++++++++++++++++!	  
999 if((exist_flag1).and.(exist_flag2)) write(6,*) 'reading chninfo.data and skip sidechain adjustment'
      	open(7,file='chninfo-n1.data',status='unknown',form='formatted') 
      do i=1,nb2
	      read(7,*) new_rxf(i),new_ryf(i),new_rzf(i)
      enddo
        CLOSE(7)	  
		
	 j = 0
          do i = 1,nb2
     			if ((i.ge.1) .and.(i.le.chnln)) then
					     j = j + 1
			       new_rx(j)=new_rxf(i)
				   new_ry(j)=new_ryf(i)
				   new_rz(j)=new_rzf(i)
				   elseif ((i.ge.(1+chnln)) .and.(i.le.(chnln+chnln))) then 
				   	     j = j + 1
			       new_rx(j)=new_rxf(i)
				   new_ry(j)=new_ryf(i)
				   new_rz(j)=new_rzf(i)
   				   elseif ((i.ge.(1+chnln*2)) .and.(i.le.(chnln+chnln*2))) then 
				   	     j = j + 1
			       new_rx(j)=new_rxf(i)
				   new_ry(j)=new_ryf(i)
				   new_rz(j)=new_rzf(i)
   				   elseif ((i.ge.(1+chnln*3)) .and.(i.le.(chnln+chnln*3))) then 
				   	     j = j + 1
			       new_rx(j)=new_rxf(i)
				   new_ry(j)=new_ryf(i)
				   new_rz(j)=new_rzf(i)
               endif
             enddo

      	open(7,file='chninfo-n2.data',status='unknown',form='formatted') 
      do i=1,nb4
	      read(7,*) new_rxf2(i),new_ryf2(i),new_rzf2(i)
      enddo
        CLOSE(7)
	 j = 0
          do i = 1,nb4
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
      do i=1, nb2-1
      do j=i+1, nb2
		dislocx(j) = new_rx(j) - new_rx(i)
		dislocy(j) = new_ry(j) - new_ry(i)
		dislocz(j) = new_rz(j) - new_rz(i)
      enddo
      enddo
      do i=1,nb2
		xtmp(i)=new_rx(i)
		ytmp(i)=new_ry(i)
		ztmp(i)=new_rz(i)
      enddo
      do i=1, nb4-1
      do j=i+1, nb4
		dislocx2(j) = new_rx2(j) - new_rx2(i)
		dislocy2(j) = new_ry2(j) - new_ry2(i)
		dislocz2(j) = new_rz2(j) - new_rz2(i)
      enddo
      enddo
      do i=1,nb4
		xtmp(nop2+i)=new_rx2(i)
		ytmp(nop2+i)=new_ry2(i)
		ztmp(nop2+i)=new_rz2(i)
      enddo

      	open(7,file='inputs/beadwell_ha55a.data',status='unknown')
      	do i=1,400
      		read(7,711) iiii,jjjj,bdtemp,wltemp
711  		format(2(2x,i2,2x),2(f6.3,2x))
       		bds(iiii,jjjj)=bdtemp
       		wel(iiii,jjjj)=wltemp
      	end do
      close(7)

      if(bds(29,29)/2 .gt. 2.000) then
		sig_max=2.000+bds(29,29)/2
      else
		sig_max=4.0D0
      endif
      write(6,*)'Maximum Sigma',sig_max
      do i=1,29
        do j=1,29
			sig_ij=bds(i,j)
	       	if (sig_ij .gt. sig_max) then
                  		sig_max = sig_ij
      write(6,*)'Maximum Sigma',sig_max
              	endif
      enddo
      enddo
		
	!this selects a random position, build the peptide in that location, and checks for overlaps
	!if there is an overlap it repeats
      do i=1,nc
		!write(6,*)'Peptide 1',i
993		xtmp((i-1)*nb2+1)=(drandm(0))*boxl
		xtmp((i-1)*nb2+1)=xtmp((i-1)*nb2+1)-boxl*ANINT(xtmp((i-1)*nb2+1)/boxl)
		ytmp((i-1)*nb2+1)=(drandm(0))*boxl
		ytmp((i-1)*nb2+1)=ytmp((i-1)*nb2+1)-boxl*ANINT(ytmp((i-1)*nb2+1)/boxl)
		ztmp((i-1)*nb2+1)=(drandm(0))*boxl
		ztmp((i-1)*nb2+1)=ztmp((i-1)*nb2+1)-boxl*ANINT(ztmp((i-1)*nb2+1)/boxl)
       do j=2,nb2
			xtmp((i-1)*nb2+j)=xtmp((i-1)*nb2+j-1)+dislocx(j)
			xtmp((i-1)*nb2+j)=xtmp((i-1)*nb2+j)-boxl*ANINT(xtmp((i-1)*nb2+j)/boxl)
			ytmp((i-1)*nb2+j)=ytmp((i-1)*nb2+j-1)+dislocy(j)
			ytmp((i-1)*nb2+j)=ytmp((i-1)*nb2+j)-boxl*ANINT(ytmp((i-1)*nb2+j)/boxl)
			ztmp((i-1)*nb2+j)=ztmp((i-1)*nb2+j-1)+dislocz(j)
			ztmp((i-1)*nb2+j)=ztmp((i-1)*nb2+j)-boxl*ANINT(ztmp((i-1)*nb2+j)/boxl)
       enddo
        do k=(i-1)*nb2+1,i*nb2
	      		do l=1,(i-1)*nb2
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
9993	xtmp(nop2+(i-1)*nb4+1)=(drandm(0))*boxl
		xtmp(nop2+(i-1)*nb4+1)=xtmp(nop2+(i-1)*nb4+1)-boxl*ANINT(xtmp(nop2+(i-1)*nb4+1)/boxl)
		ytmp(nop2+(i-1)*nb4+1)=(drandm(0))*boxl
		ytmp(nop2+(i-1)*nb4+1)=ytmp(nop2+(i-1)*nb4+1)-boxl*ANINT(ytmp(nop2+(i-1)*nb4+1)/boxl)
		ztmp(nop2+(i-1)*nb4+1)=(drandm(0))*boxl
		ztmp(nop2+(i-1)*nb4+1)=ztmp(nop2+(i-1)*nb4+1)-boxl*ANINT(ztmp(nop2+(i-1)*nb4+1)/boxl)
       do j=2,nb4
			xtmp(nop2+(i-1)*nb4+j)=xtmp(nop2+(i-1)*nb4+j-1)+dislocx2(j)
			xtmp(nop2+(i-1)*nb4+j)=xtmp(nop2+(i-1)*nb4+j)-boxl*ANINT(xtmp(nop2+(i-1)*nb4+j)/boxl)
			ytmp(nop2+(i-1)*nb4+j)=ytmp(nop2+(i-1)*nb4+j-1)+dislocy2(j)
			ytmp(nop2+(i-1)*nb4+j)=ytmp(nop2+(i-1)*nb4+j)-boxl*ANINT(ytmp(nop2+(i-1)*nb4+j)/boxl)
			ztmp(nop2+(i-1)*nb4+j)=ztmp(nop2+(i-1)*nb4+j-1)+dislocz2(j)
			ztmp(nop2+(i-1)*nb4+j)=ztmp(nop2+(i-1)*nb4+j)-boxl*ANINT(ztmp(nop2+(i-1)*nb4+j)/boxl)
       enddo
        do k=(i-1)*nb4+1,i*nb4
	      		do l=1,(i-1)*nb4
                       		rxij=xtmp(nop2+k)-xtmp(nop2+l)
                       		ryij=ytmp(nop2+k)-ytmp(nop2+l) 
                       		rzij=ztmp(nop2+k)-ztmp(nop2+l)
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
      open(7,file='checks/config.pdb')
      do i=1,nop2+nop4
       write(7,'(a6,i5,a3,1x,15x,3f8.3)') 'ATOM  ',i,'N', xtmp(i), ytmp(i), ztmp(i)
      enddo
      close(7)

!      open(7,file='checks/emblem.input')
!	j=0
!      do i=1,chnln
!	j=j+1
!		write(7,'(2i5,15x,3f8.3,i5)') j,j, xtmp(chnln+i), ytmp(chnln+i), ztmp(chnln+i),identity(chnln+i)
!	j=j+1
!		write(7,'(2i5,15x,3f8.3,i5)') j,j, xtmp(i), ytmp(i), ztmp(i), identity(i)
!	j=j+1
!		write(7,'(2i5,15x,3f8.3,i5)') j,j, xtmp(chnln*2+i), ytmp(chnln*2+i), ztmp(chnln*2+i), identity(chnln*2+i)
!	j=j+1
!		write(7,'(2i5,15x,3f8.3,i5)') j,j, xtmp(chnln*3+i), ytmp(chnln*3+i), ztmp(chnln*3+i), identity(chnln*3+i)
!      enddo
!         close(7)
   	
      do i=1,nop2+nop4
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
      open(7,file='checks/configcheck.out')
      write(7,*) coll,t
      write(7,*)xtmp/boxl,ytmp/boxl,ztmp/boxl
      close(7)	
	
	!this writes the actual config input
      open(7,file='results/run0000.config',status='unknown',form='unformatted')
        write(7) coll,t,xtmp/boxl,ytmp/boxl,ztmp/boxl
      close(7)
	
	!read in masses and adjust for bead identity
      open(7,file='parameters/mass.data',status='old')
      do i=1,23
       read(7,'(4x,i2,2x,f8.3)') iiii,bmass_temp
        bmass(iiii)=bmass_temp
      end do
	
	bmass(3)=bmass(20)
	
      do i=1,4
		bmass(i+4)=bmass(i)
      end do

      close(7)

      do i=1,nop
		bm(i)=bmass(identity(i))
      end do
      do i=1,nop3
		 bm2(i)=bmass(identity2(i))
      end do	  

	!write to make sure masses are correct
      open(7,file='checks/masses.out')
       do i=1,nop
        write(7,*) i,identity(i),bm(i)
      enddo
       do i=1,nop3
        write(7,*) nop+i,identity2(i),bm2(i)
      enddo
       close(7)

      write(6,*)'Desired reduced temperature',setemp/12

	!assign random velocities
      do k=1,nop
       if(identity(k).ne.9) then
         	x1=drandm(0)*0.75d0+0.2d0
         	x2=drandm(0)*0.75d0+0.2d0
         	y1=drandm(0)*0.75d0+0.2d0
         	y2=drandm(0)*0.75d0+0.2d0
         	sv(4,k)=dsqrt(-2.d0*log(x1))*cos(2.d0*pi*y1)
         	sv(5,k)=dsqrt(-2.d0*log(x1))*sin(2.d0*pi*y1)
         	sv(6,k)=dsqrt(-2.d0*log(x2))*cos(2.d0*pi*y2)
       endif
      end do
      do k=1,nop3
       if(identity2(k).ne.9) then
         	x1=drandm(0)*0.75d0+0.2d0
         	x2=drandm(0)*0.75d0+0.2d0
         	y1=drandm(0)*0.75d0+0.2d0
         	y2=drandm(0)*0.75d0+0.2d0
         	sv(4,nop+k)=dsqrt(-2.d0*log(x1))*cos(2.d0*pi*y1)
         	sv(5,nop+k)=dsqrt(-2.d0*log(x1))*sin(2.d0*pi*y1)
         	sv(6,nop+k)=dsqrt(-2.d0*log(x2))*cos(2.d0*pi*y2)
       endif
      end do

	!scale velocities so that linear momentum = 0
      do k=1,nop
       if(identity(k).ne.9) then
         	sumx=sumx+sv(4,k)
         	sumy=sumy+sv(5,k)
         	sumz=sumz+sv(6,k)
      		endif
      end do
      do k=1,nop3
       if(identity2(k).ne.9) then
         	sumx=sumx+sv(4,nop+k)
         	sumy=sumy+sv(5,nop+k)
         	sumz=sumz+sv(6,nop+k)
      		endif
      end do	
	
      do k=1,nop
       if(identity(k).ne.9) then
         	sv(4,k)=sv(4,k)-sumx/real(nop2+nop4)
         	sv(5,k)=sv(5,k)-sumy/real(nop2+nop4)
         	sv(6,k)=sv(6,k)-sumz/real(nop2+nop4)
     	endif
      end do
      do k=1,nop3
       if(identity2(k).ne.9) then
         	sv(4,nop+k)=sv(4,nop+k)-sumx/real(nop2+nop4)
         	sv(5,nop+k)=sv(5,nop+k)-sumy/real(nop2+nop4)
         	sv(6,nop+k)=sv(6,nop+k)-sumz/real(nop2+nop4)
     	endif
      end do	
	sumvel=0.d0
	
	!check summed velocities
      open(7,file='checks/sumvelcheck.out')
      do k=1,nop
       if(identity(k).ne.9) then
			sumvel=sumvel+(sv(4,k)**2+sv(5,k)**2+sv(6,k)**2)/bm(k)
       endif
      enddo  
      do k=1,nop3
       if(identity2(k).ne.9) then
			sumvel=sumvel+(sv(4,nop+k)**2+sv(5,nop+k)**2+sv(6,nop+k)**2)/bm2(k)
       endif
      enddo 
	!here sv is actually momentum not velocity
	 
	tred=sumvel/3.d0/dble(nop2+nop4)
	const=dsqrt(setemp/tred)

       do k=1,nop
       if(identity(k).ne.9) then
     		sv(4,k)=const*sv(4,k)/bm(k)
     		sv(5,k)=const*sv(5,k)/bm(k)
     		sv(6,k)=const*sv(6,k)/bm(k)
        endif
      enddo
       do k=1,nop3
       if(identity2(k).ne.9) then
     		sv(4,nop+k)=const*sv(4,nop+k)/bm2(k)
     		sv(5,nop+k)=const*sv(5,nop+k)/bm2(k)
     		sv(6,nop+k)=const*sv(6,nop+k)/bm2(k)
        endif
      enddo
  	
	sumvel=0.d0
      
      do k=1,nop
        if(identity(k).ne.9) then
			sumvel=sumvel+bm(k)*(sv(4,k)**2+sv(5,k)**2+sv(6,k)**2)
			write(7,*)k,sumvel
      	endif
      enddo
      do k=1,nop3
        if(identity2(k).ne.9) then
			sumvel=sumvel+bm2(k)*(sv(4,nop+k)**2+sv(5,nop+k)**2+sv(6,nop+k)**2)
			write(7,*)nop+k,sumvel
      	endif
      enddo
      close(7)
	
	tred=sumvel/3.d0/dble(nop2+nop4)
      write(6,*)'New reduced temperature',tred/12.0

	!this is the velocity input to the simulations
      open(7,file='results/run0000.lastvel',status='unknown',form='unformatted') 
	l=0
      do k=1,nop   
         if(identity(k).ne.9) then
         	l=l+1
         	xtmp(l)=sv(4,k)
         	ytmp(l)=sv(5,k)
         	ztmp(l)=sv(6,k)
        endif
      enddo
      do k=1,nop3   
         if(identity2(k).ne.9) then
         	l=l+1
         	xtmp(l)=sv(4,nop+k)
         	ytmp(l)=sv(5,nop+k)
         	ztmp(l)=sv(6,nop+k)
        endif
      enddo
      write(7) coll,xtmp,ytmp,ztmp
      close(7)

	!these are the two blank files needed to start
      open(7,file='results/run0000.bptnr',status='unknown',form='unformatted')
      close(7)
      open(7,file='results/run0000.energy',status='unknown',form='unformatted')
      close(7)

	!check velocities
      open(7,file='checks/velocitycheck.out')
	l=0
       do k=1,nop  
       if(identity(k).ne.9) then
         	l=l+1
         	xtmp(l)=sv(4,k)
         	ytmp(l)=sv(5,k)
         	ztmp(l)=sv(6,k)
        write(7,*) k,identity(k),xtmp(k),ytmp(k),ztmp(k)
    		endif
      enddo
       do k=1,nop3   
       if(identity2(k).ne.9) then
         	l=l+1
         	xtmp(l)=sv(4,nop+k)
         	ytmp(l)=sv(5,nop+k)
         	ztmp(l)=sv(6,nop+k)
        write(7,*) nop+k,identity2(k),xtmp(nop+k),ytmp(nop+k),ztmp(nop+k)
    		endif
      enddo
      close(7)

	!check HB parnters
      open(7,file='checks/hbcheck.out')
      do i=1,nop2
		write(7,*)i,bptnr(i)
      enddo
      do i=1,nop4
		write(7,*)nop+i,bptnr2(i)
      enddo	  
      close(7)

      stop
      end
