        program removepbc
		
       implicit none


        integer nb,nc,chnln,nop

        parameter(nb=55)
        parameter(chnln=11)
        parameter(nc=40)
        parameter(nop=nb*nc)
		
        character*1 col(nop*8),newcol(nop*8),xcol(nb*8),newcol2(nop*8),ncol(nop)
        character*3 amino_name(nop*8),newamino_name(nop*8),xamino_name(nb*8),newamino_name2(nop*8),namino_name(nop)
        character*4 elementtype(nop*8),newelementtype(nop*8),xelementtype(nb*8),newelementtype2(nop*8),nelementtype(nop),atom(nop*8),newatom(nop*8),xatom(nb*8),newatom2(nop*8),natom(nop)
        integer check,a,b,c,i,j,k,l,m,n,p,label(nop*8),newlabel(nop*8),xlabel(nb*8),newlabel2(nop*8),nlabel(nop),res(nop*8),newres(nop*8),xres(nb*8),newres2(nop*8),nres(nop),chnnum(nop)
		integer i1,i2,i3,i4,i5,i6,i7
		integer newcount(nc),newcountx(nc),newcounty(nc),newcountz(nc),count(nc),count1,countx(nc),county(nc),countz(nc),target1(nc,3),target2(nc,3),list(nc*2),fside1(chnln)
		real*8 rx(nop*8),ry(nop*8),rz(nop*8),newrx(nop*8),newry(nop*8),newrz(nop*8),xrx(nb*8),xry(nb*8),xrz(nb*8),boxl,newrx2(nop*8),newry2(nop*8),newrz2(nop*8),nrx(nop),nry(nop),nrz(nop)
		real*8 increx,increy,increz
		!real*8 wcx,wcy,wcz,d1,d2,c1x(nb),c1y(nb),c1z(nb)
		real*8 avgnewrx,avgnewry,avgnewrz,avgc1x,avgc1y,avgc1z
		real*8 bondx(chnln*nc,5),bondy(chnln*nc,5),bondz(chnln*nc,5),bond(chnln*nc,5)
		real*8 newbondx(chnln*nc,5),newbondy(chnln*nc,5),newbondz(chnln*nc,5),newbond(chnln*nc,5)
		!real*8 xp,xn,yp,yn,zp,zn
		real*8 com1,com2,com3,com4
        real*8 factorx,factory,factorz
		integer sumcount,sumnewcount
		     logical jx(nc),jy(nc),jz(nc) 
			 
			 
		character*64 input1,output1,input2  !,suffix
        character*4 fname_digits
      common fname_digits		
	integer izero,iflag
	character*64 filename
	logical exist_flag
	character*1 zero
	character*1 char
	integer ichar, len
	data zero/'0'/
	
	izero = ichar('0')
	iflag = 0
	filename = '../run'//'0000.config'
	fname_digits = '0000'
	inquire(file = filename, exist = exist_flag)

	7       format(A4,3X,I4,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3)
         !real*8 dx,dy,dz,dr
		 
		!read(5,*) suffix			
			
			boxl=149.0d0

	    input2='../../parameters/firstside1.data'
			open(unit=7,file=input2,action='read',status='unknown')
      	read(7,*) fside1
      	close(unit=7)
      	l=1	
      	do k=1,chnln
         	if (fside1(k) .ne. 0) then
	    		fside1(k)=3*chnln+l
	    		!write(6,*) k, fside1(k)
	    		l=l+1
				endif
      	enddo
	
	
	do while (exist_flag)		
		iflag = iflag + 1

		fname_digits = char(iflag/1000+izero)//char(mod(iflag,1000)/100+izero)//char(mod(iflag,100)/10+izero)//char(mod(iflag,10)+izero)
	    input1='../run'//fname_digits//'.pdb'
	    output1=''//fname_digits//''
	   	inquire( file = input1, exist = exist_flag)
		write(6,*) 'iflag: ',iflag
		
		
				 do i=1,nc
					jx(i)=.false.
					jy(i)=.false.
					jz(i)=.false.
				  enddo
				  
      !suffix='009.pdb'			   	  
	  l=0
	  
		open(9,file=input1,action='read',form='formatted')
					   do i=1,nc
		       count(i) = 0
		       countx(i) = 0
		       county(i) = 0
		       countz(i) = 0
			   enddo
			   
		       do k=1,nc
			  ! m =(k-1)*(nb2+chnln)+1
			   m=(k-1)*nb +1
				j =1		  
				p=0
 					do while (j .le. chnln)

                  read(9,7) atom(m),label(m),elementtype(m),amino_name(m),col(m),res(m), rx(m), ry(m), rz(m)							
				  read(9,7) atom(m+1),label(m+1),elementtype(m+1),amino_name(m+1),col(m+1),res(m+1), rx(m+1), ry(m+1), rz(m+1)		
				  read(9,7) atom(m+2),label(m+2),elementtype(m+2),amino_name(m+2),col(m+2),res(m+2), rx(m+2), ry(m+2), rz(m+2)			
				  read(9,7) atom(m+3),label(m+3),elementtype(m+3),amino_name(m+3),col(m+3),res(m+3), rx(m+3), ry(m+3), rz(m+3)	

				  rx(m+3)=rx(m+3)-boxl*dnint(rx(m+3)/boxl)		
					ry(m+3)=ry(m+3)-boxl*dnint(ry(m+3)/boxl)
					rz(m+3)=rz(m+3)-boxl*dnint(rz(m+3)/boxl)	
					
				  !1:N-CA, 2:CA-CO 3:CO=O 4 :R-CA 5: CO(n-1)-N
				  p = (k-1)*chnln+j   ! p is the residue label
				 bondx(p,1)=rx(m)-rx(m+1)
				 bondy(p,1)=ry(m)-ry(m+1)
				 bondz(p,1)=rz(m)-rz(m+1)
				 bondx(p,2)=rx(m+1)-rx(m+2)
				 bondy(p,2)=ry(m+1)-ry(m+2)
				 bondz(p,2)=rz(m+1)-rz(m+2)
				  bondx(p,3)=rx(m+2)-rx(m+3)
				 bondy(p,3)=ry(m+2)-ry(m+3)
				 bondz(p,3)=rz(m+2)-rz(m+3)
				 bond(p,1)=dsqrt(bondx(p,1)**2+bondy(p,1)**2+bondz(p,1)**2)
				 bond(p,2)=dsqrt(bondx(p,2)**2+bondy(p,2)**2+bondz(p,2)**2)
				 bond(p,3)=dsqrt(bondx(p,3)**2+bondy(p,3)**2+bondz(p,3)**2)
				 
				 !if(k.eq.1) then
				!	write(6,*)'bond: ',k,j,bond(p,3)
			    ! endif
				 
				  if (j.ne.1) then
 				 bondx(p,5)=rx(b)-rx(m)
				 bondy(p,5)=ry(b)-ry(m)
				 bondz(p,5)=rz(b)-rz(m)
				 bond(p,5)=dsqrt(bondx(p,5)**2+bondy(p,5)**2+bondz(p,5)**2)		
				 endif
		
			    b=m+2  ! restore the CO label of the previous residue 
				 
				 
	    if (fside1(j).ne.0) then
			  read(9,7) atom(m+4),label(m+4),elementtype(m+4),amino_name(m+4),col(m+4),res(m+4), rx(m+4), ry(m+4), rz(m+4)				

 				 bondx(p,4)=rx(m+4)-rx(m+1)
				 bondy(p,4)=ry(m+4)-ry(m+1)
				 bondz(p,4)=rz(m+4)-rz(m+1)
				 bond(p,4)=dsqrt(bondx(p,4)**2+bondy(p,4)**2+bondz(p,4)**2)		
				
				 m = m + 5		
		else 			
                 m = m + 4			
		endif	
		          l=0
				   do a=1,5
				   
				   if(bond(p,a).ge.boxl/2) then
				    l=l+1
					list(l)=k
			       
				   !write(6,*)'# bond,chain,residue,bty,blength: ', l,k,p,a,bond(p,a)
				   
	               if(abs(bondx(p,a)).ge.boxl*0.9) then
					 jx(k)=.true.
					 !target1(l,1)=m
					 !target2(l,1)=n
					 countx(k)=1                ! # bead of a chain penetrating wall
					 endif
				   if(abs(bondy(p,a)).ge.boxl*0.9) then
					 jy(k)=.true.
					 !target1(l,2)=m
					 !target2(l,2)=n
 					 county(k)=1
					  endif
                   if(abs(bondz(p,a)).ge.boxl*0.9) then	
					 jz(k)=.true.	
					 !target1(l,3)=m
					 !target2(l,3)=n					 
					 countz(k)=1					 
				    endif				   
				   
				    endif
					
                 enddo
				   
				     j = j + 1
					 
					 enddo  		 
				count(k)=countx(k)+county(k)+countz(k)			
                   ! write(6,*)'# of chain penetrating boudnary: ',count(k),countx(k),county(k),countz(k)				
			   enddo		


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
   write(6,*)'start moving the system as a whole part'

		
		do factorx=0.0,boxl,10.0d0
		do factory=0.0,boxl,10.0d0		
		do factorz=0.0,boxl,10.0d0		
		       sumnewcount = 0
			   do i=1,nc
		       newcount(i) = 0
		       newcountx(i) = 0
		       newcounty(i) = 0
		       newcountz(i) = 0
			   enddo
		do i=1,nop
		
	 newatom(i)=atom(i)
 newlabel(i)=label(i)
 newelementtype(i)=elementtype(i)
 newamino_name(i)=amino_name(i)
 newcol(i)=col(i)
 newres(i)=res(i)
 newrx(i)=rx(i)
 newry(i)=ry(i)
 newrz(i)=rz(i)	
 
		i1=nop+i
		i2=nop*2+i
		i3=nop*3+i
		i4=nop*4+i
		i5=nop*5+i
		i6=nop*6+i
		i7=nop*7+i
		
 newatom(i1)=atom(i)
 newlabel(i1)=label(i)
 newelementtype(i1)=elementtype(i)
 newamino_name(i1)=amino_name(i)
 newcol(i1)=col(i)
 newres(i1)=res(i)
 newrx(i1)=rx(i)+boxl
 newry(i1)=ry(i)
 newrz(i1)=rz(i)	

  newatom(i2)=atom(i)
 newlabel(i2)=label(i)
 newelementtype(i2)=elementtype(i)
 newamino_name(i2)=amino_name(i)
 newcol(i2)=col(i)
 newres(i2)=res(i)
 newrx(i2)=rx(i)
 newry(i2)=ry(i)+boxl
 newrz(i2)=rz(i)
 
  newatom(i3)=atom(i)
 newlabel(i3)=label(i)
 newelementtype(i3)=elementtype(i)
 newamino_name(i3)=amino_name(i)
 newcol(i3)=col(i)
 newres(i3)=res(i)
 newrx(i3)=rx(i)
 newry(i3)=ry(i)
 newrz(i3)=rz(i)+boxl
 
  newatom(i4)=atom(i)
 newlabel(i4)=label(i)
 newelementtype(i4)=elementtype(i)
 newamino_name(i4)=amino_name(i)
 newcol(i4)=col(i)
 newres(i4)=res(i)
 newrx(i4)=rx(i)+boxl
 newry(i4)=ry(i)+boxl
 newrz(i4)=rz(i)	

  newatom(i5)=atom(i)
 newlabel(i5)=label(i)
 newelementtype(i5)=elementtype(i)
 newamino_name(i5)=amino_name(i)
 newcol(i5)=col(i)
 newres(i5)=res(i)
 newrx(i5)=rx(i)
 newry(i5)=ry(i)+boxl
 newrz(i5)=rz(i)+boxl
 
  newatom(i6)=atom(i)
 newlabel(i6)=label(i)
 newelementtype(i6)=elementtype(i)
 newamino_name(i6)=amino_name(i)
 newcol(i6)=col(i)
 newres(i6)=res(i)
 newrx(i6)=rx(i)+boxl
 newry(i6)=ry(i)
 newrz(i6)=rz(i)+boxl
 
 
  newatom(i7)=atom(i)
 newlabel(i7)=label(i)
 newelementtype(i7)=elementtype(i)
 newamino_name(i7)=amino_name(i)
 newcol(i7)=col(i)
 newres(i7)=res(i)
 newrx(i7)=rx(i)+boxl
 newry(i7)=ry(i)+boxl
 newrz(i7)=rz(i)+boxl	
 
		enddo
		
		
		do i=1,8*nop
!#ifdef  bias
		if((newrx(i).le.(boxl/2+factorx)).and.(newrx(i).ge.(-boxl/2+factorx)).and.(newry(i).le.(boxl/2+factory)).and.(newry(i).ge.(-boxl/2+factory)) .and.(newrz(i).le.(boxl/2+factorz)).and.(newrz(i).ge.(-boxl/2+factorz))) then
!#else
!		if((newrx(i).le.(boxl)).and.(newrx(i).ge.(0.0)).and.(newry(i).le.(boxl)).and.(newry(i).ge.(0.0)) .and.(newrz(i).le.(boxl)).and.(newrz(i).ge.(0.0))) then
!#endif
		!count1=count1+1
			!write(6,*)'i: ',i
		    if(int(MOD(i,nop)).ne.0) then 
			m=label(int(MOD(i,nop)))
			else
			!write(6,*)'the last bead!'
			m=label(nop)
			endif

          newatom2(m)=newatom(i)
		  newlabel2(m)=newlabel(i)
		  newelementtype2(m)=newelementtype(i)
		  newamino_name2(m)=newamino_name(i)
		  newcol2(m)=newcol(i)
		  newres2(m)=newres(i)
		  newrx2(m)=newrx(i)
	 	  newry2(m)=newry(i)
		  newrz2(m)=newrz(i)			
		
		endif
		
		enddo
		!write(6,*)'in 2nd part: nop and count: ',nop,count1

		do i=1,nop
		  newrx2(i)=newrx2(i)-factorx !+boxl/2
	 	  newry2(i)=newry2(i)-factory !+boxl/2
		  newrz2(i)=newrz2(i)-factorz !+boxl/2
		enddo	
			
		close(9)		
		
		
		! last procedure: check if there is any penetrating boundary happens!  
		  
		   do k=1,nc
			   m=(k-1)*nb +1
				j =1		  
				p=0
 					do while (j .le. chnln)
					
				  !1:N-CA, 2:CA-CO 3:CO=O 4 :R-CA 5: CO(n-1)-N
				  p = (k-1)*chnln+j   ! p is the residue label
				 newbondx(p,1)=newrx2(m)-newrx2(m+1)
				 newbondy(p,1)=newry2(m)-newry2(m+1)
				 newbondz(p,1)=newrz2(m)-newrz2(m+1)
				 newbondx(p,2)=newrx2(m+1)-newrx2(m+2)
				 newbondy(p,2)=newry2(m+1)-newry2(m+2)
				 newbondz(p,2)=newrz2(m+1)-newrz2(m+2)
				  newbondx(p,3)=newrx2(m+2)-newrx2(m+3)
				 newbondy(p,3)=newry2(m+2)-newry2(m+3)
				 newbondz(p,3)=newrz2(m+2)-newrz2(m+3)
				 newbond(p,1)=dsqrt(newbondx(p,1)**2+newbondy(p,1)**2+newbondz(p,1)**2)
				 newbond(p,2)=dsqrt(newbondx(p,2)**2+newbondy(p,2)**2+newbondz(p,2)**2)
				 newbond(p,3)=dsqrt(newbondx(p,3)**2+newbondy(p,3)**2+newbondz(p,3)**2)
				 
				  if (j.ne.1) then
 				 newbondx(p,5)=newrx2(b)-newrx2(m)
				 newbondy(p,5)=newry2(b)-newry2(m)
				 newbondz(p,5)=newrz2(b)-newrz2(m)
				 newbond(p,5)=dsqrt(newbondx(p,5)**2+newbondy(p,5)**2+newbondz(p,5)**2)		
				   endif
		
			    b=m+2  ! restore the CO label of the previous residue 
				 		 
	    if (fside1(j).ne.0) then
		
 				 newbondx(p,4)=newrx2(m+4)-newrx2(m+1)
				 newbondy(p,4)=newry2(m+4)-newry2(m+1)
				 newbondz(p,4)=newrz2(m+4)-newrz2(m+1)
				 newbond(p,4)=dsqrt(newbondx(p,4)**2+newbondy(p,4)**2+newbondz(p,4)**2)		
				
				 m = m + 5		
		else 			
                 m = m + 4			
		endif	
		          l=0
				  
				   do a=1,5
				   
				   !if(newbond(p,a).ge.boxl/2) then			       
				   !write(6,*)'# newbond,chain,residue,bty,blength: ', l,k,p,a,newbond(p,a)
				   
	               if(abs(newbondx(p,a)).ge.boxl*0.9) then
					 newcountx(k)=1                ! # bead of a chain penetrating wall
					 endif
				   if(abs(newbondy(p,a)).ge.boxl*0.9) then
 					 newcounty(k)=1
					  endif
                   if(abs(newbondz(p,a)).ge.boxl*0.9) then			 
					 newcountz(k)=1					 
				    endif				   				   
				    !endif					
                 enddo
				 	   
				     j = j + 1
					 
					 enddo  		 
					 newcount(k)=newcountx(k)+newcounty(k)+newcountz(k)
				sumnewcount=sumnewcount+newcount(k)		
					
			   enddo		
		
					   if(sumnewcount.eq.0) then	   

	 	      open(9,file=output1,action='write',form='formatted')	
		       do i=1,nop
		 write(9,7) newatom2(i),newlabel2(i),newelementtype2(i),newamino_name2(i),newcol2(i),newres2(i), newrx2(i), newry2(i), newrz2(i)	
		        enddo	
		    close(9)
                      goto 100					  
			        endif	
		
		enddo
		enddo
		enddo
		
100         if(sumnewcount.eq.0) then
             write(6,*)'NO chain penetrates boundary !'	
			 
			    goto 200
			 
     		else
			  write(6,*)'chain penetrates boundary !', sumnewcount	
			endif   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

increx=0.0
increy=0.0
increz=0.0
			do c=1,nc

				if(count(c).eq.0)  then
				do j=1,nb
				      i=(c-1)*nb+j
					newatom(i)=atom(i)
					newlabel(i)=label(i)
					newelementtype(i)=elementtype(i)
					newamino_name(i)=amino_name(i)
					newcol(i)=col(i)
					newres(i)=res(i)
					newrx(i)=rx(i)
					newry(i)=ry(i)
					newrz(i)=rz(i)	
                  enddo					
				endif
		     enddo

	 l = 0
		
		do j=1,nc
		      
			  if(count(j).eq.1) then		  
	
		     do c=1,nb
		     i=(j-1)*nb+c
		       i1=nb+c

			    xatom(c)=atom(i)
                xlabel(c)=label(i)
                xelementtype(c)=elementtype(i)
                xamino_name(c)=amino_name(i)
                xcol(c)=col(i)
                xres(c)=res(i)
                xrx(c)=rx(i)
                xry(c)=ry(i)
                xrz(c)=rz(i)	
			   
			  if(jx(j)) then
			  increx=boxl/2
			  increy=0.0
			  increz=0.0			  
 xatom(i1)=atom(i)
 xlabel(i1)=label(i)
 xelementtype(i1)=elementtype(i)
 xamino_name(i1)=amino_name(i)
 xcol(i1)=col(i)
 xres(i1)=res(i)
 xrx(i1)=rx(i)+boxl
 xry(i1)=ry(i)
 xrz(i1)=rz(i)		 
			  elseif(jy(j))then
	          increy=boxl/2
			  increx=0.0
			  increz=0.0
 xatom(i1)=atom(i)
 xlabel(i1)=label(i)
 xelementtype(i1)=elementtype(i)
 xamino_name(i1)=amino_name(i)
 xcol(i1)=col(i)
 xres(i1)=res(i)
 xrx(i1)=rx(i)
 xry(i1)=ry(i)+boxl
 xrz(i1)=rz(i)
 
			  elseif(jz(j))then
         increz=boxl/2
			  increx=0.0
			  increy=0.0
		 
 xatom(i1)=atom(i)
 xlabel(i1)=label(i)
 xelementtype(i1)=elementtype(i)
 xamino_name(i1)=amino_name(i)
 xcol(i1)=col(i)
 xres(i1)=res(i)
 xrx(i1)=rx(i)
 xry(i1)=ry(i)
 xrz(i1)=rz(i)+boxl
              endif 
			  
              enddo			  

          check=0
		  
	do i=1,2*nb
			
		if((xrx(i).le.(boxl/2+increx)).and.(xrx(i).ge.(-boxl/2+increx)).and.(xry(i).le.(boxl/2+increy)).and.(xry(i).ge.(-boxl/2+increy)) .and.(xrz(i).le.(boxl/2+increz)).and.(xrz(i).ge.(-boxl/2+increz))) then
			
			check=check+1
			!if(j.lt.20) then
			!	write(6,*)'check chain: ',j, check
			!endif
		    if(int(MOD(i,nb)).ne.0) then 
			   m=xlabel(int(MOD(i,nb)))
			else
			  ! write(6,*)'1st case, the last bead!',j
			   m=xlabel(nb)
			endif		
			
          newatom(m)=xatom(i)
		  newlabel(m)=xlabel(i)
		  newelementtype(m)=xelementtype(i)
		  newamino_name(m)=xamino_name(i)
		  newcol(m)=xcol(i)
		  newres(m)=xres(i)
		  newrx(m)=xrx(i)
	 	  newry(m)=xry(i)
		  newrz(m)=xrz(i)		
	
			endif
		    enddo	  		
			
 !        write(6,*)'check must be equal to: ',nb,check

	! finish the first case
	
			  elseif(count(j).eq.2) then		  
             check=0
			  
		      do c=1,nb
		      i=(j-1)*nb+c
		      i1=nb+c
		      i2=nb*2+c
		      i3=nb*3+c
	
			    xatom(c)=atom(i)
                xlabel(c)=label(i)
                xelementtype(c)=elementtype(i)
                xamino_name(c)=amino_name(i)
                xcol(c)=col(i)
                xres(c)=res(i)
                xrx(c)=rx(i)
                xry(c)=ry(i)
                xrz(c)=rz(i)	

	
			  if(jx(j).and.jy(j)) then
			  increx=boxl/2
			  increy=boxl/2
               increz=0.0			  
			  
 xatom(i1)=atom(i)
 xlabel(i1)=label(i)
 xelementtype(i1)=elementtype(i)
 xamino_name(i1)=amino_name(i)
 xcol(i1)=col(i)
 xres(i1)=res(i)
 xrx(i1)=rx(i)+boxl
 xry(i1)=ry(i)
 xrz(i1)=rz(i)	
 
 xatom(i2)=atom(i)
 xlabel(i2)=label(i)
 xelementtype(i2)=elementtype(i)
 xamino_name(i2)=amino_name(i)
 xcol(i2)=col(i)
 xres(i2)=res(i)
 xrx(i2)=rx(i)
 xry(i2)=ry(i)+boxl
 xrz(i2)=rz(i)	
 
 
 xatom(i3)=atom(i)
 xlabel(i3)=label(i)
 xelementtype(i3)=elementtype(i)
 xamino_name(i3)=amino_name(i)
 xcol(i3)=col(i)
 xres(i3)=res(i)
 xrx(i3)=rx(i)+boxl
 xry(i3)=ry(i)+boxl
 xrz(i3)=rz(i)
 
			  elseif(jy(j).and.jz(j))then
			  increx=0.0
			  increy=boxl/2
			  increz=boxl/2	
 xatom(i1)=atom(i)
 xlabel(i1)=label(i)
 xelementtype(i1)=elementtype(i)
 xamino_name(i1)=amino_name(i)
 xcol(i1)=col(i)
 xres(i1)=res(i)
 xrx(i1)=rx(i)
 xry(i1)=ry(i)+boxl
 xrz(i1)=rz(i)	
 			  
 xatom(i2)=atom(i)
 xlabel(i2)=label(i)
 xelementtype(i2)=elementtype(i)
 xamino_name(i2)=amino_name(i)
 xcol(i2)=col(i)
 xres(i2)=res(i)
 xrx(i2)=rx(i)
 xry(i2)=ry(i)
 xrz(i2)=rz(i)+boxl			  
	
 xatom(i3)=atom(i)
 xlabel(i3)=label(i)
 xelementtype(i3)=elementtype(i)
 xamino_name(i3)=amino_name(i)
 xcol(i3)=col(i)
 xres(i3)=res(i)
 xrx(i3)=rx(i)
 xry(i3)=ry(i)+boxl
 xrz(i3)=rz(i)+boxl
	
			  elseif(jx(j).and.jz(j))then
			  increx=boxl/2
			  increz=boxl/2
             increy=0.0			  
 xatom(i1)=atom(i)
 xlabel(i1)=label(i)
 xelementtype(i1)=elementtype(i)
 xamino_name(i1)=amino_name(i)
 xcol(i1)=col(i)
 xres(i1)=res(i)
 xrx(i1)=rx(i)+boxl
 xry(i1)=ry(i)
 xrz(i1)=rz(i)	

 xatom(i2)=atom(i)
 xlabel(i2)=label(i)
 xelementtype(i2)=elementtype(i)
 xamino_name(i2)=amino_name(i)
 xcol(i2)=col(i)
 xres(i2)=res(i)
 xrx(i2)=rx(i)
 xry(i2)=ry(i)
 xrz(i2)=rz(i)+boxl
 
 
 xatom(i3)=atom(i)
 xlabel(i3)=label(i)
 xelementtype(i3)=elementtype(i)
 xamino_name(i3)=amino_name(i)
 xcol(i3)=col(i)
 xres(i3)=res(i)
 xrx(i3)=rx(i)+boxl
 xry(i3)=ry(i)
 xrz(i3)=rz(i)+boxl
 
			  endif
              
			  enddo					  	  
			  
	do i=1,4*nb
		if((xrx(i).le.(boxl/2+increx)).and.(xrx(i).ge.(-boxl/2+increx)).and.(xry(i).le.(boxl/2+increy)).and.(xry(i).ge.(-boxl/2+increy)) .and.(xrz(i).le.(boxl/2+increz)).and.(xrz(i).ge.(-boxl/2+increz))) then
		    check=check+1
		    if(int(MOD(i,nb)).ne.0) then 
			   m=xlabel(int(MOD(i,nb)))
			else
			  ! write(6,*)'2nd case, the last bead!',j
			   m=xlabel(nb)
			endif
			
          newatom(m)=xatom(i)
		  newlabel(m)=xlabel(i)
		  newelementtype(m)=xelementtype(i)
		  newamino_name(m)=xamino_name(i)
		  newcol(m)=xcol(i)
		  newres(m)=xres(i)
		  newrx(m)=xrx(i)
	 	  newry(m)=xry(i)
		  newrz(m)=xrz(i)	
				 		
			endif
		    enddo	  		   
! write(6,*)'check must be equal to: ',nb,check
			
			!2nd case finish
			  elseif(count(j).eq.3) then		  
             check=0
			  
		      do c=1,nb
		      i=(j-1)*nb+c
			  
			 i1=nb+c
			 i2=nb*2+c
			 i3=nb*3+c
			 i4=nb*4+c
			 i5=nb*5+c
			 i6=nb*6+c
			 i7=nb*7+c			  

			    xatom(c)=atom(i)
                xlabel(c)=label(i)
                xelementtype(c)=elementtype(i)
                xamino_name(c)=amino_name(i)
                xcol(c)=col(i)
                xres(c)=res(i)
                xrx(c)=rx(i)
                xry(c)=ry(i)
                xrz(c)=rz(i)		
				
			  increx=boxl/2
			  increy=boxl/2
              increz=boxl/2	  
			  
 xatom(i1)=atom(i)
 xlabel(i1)=label(i)
 xelementtype(i1)=elementtype(i)
 xamino_name(i1)=amino_name(i)
 xcol(i1)=col(i)
 xres(i1)=res(i)
 xrx(i1)=rx(i)+boxl
 xry(i1)=ry(i)
 xrz(i1)=rz(i)	
 
 xatom(i2)=atom(i)
 xlabel(i2)=label(i)
 xelementtype(i2)=elementtype(i)
 xamino_name(i2)=amino_name(i)
 xcol(i2)=col(i)
 xres(i2)=res(i)
 xrx(i2)=rx(i)
 xry(i2)=ry(i)+boxl
 xrz(i2)=rz(i)	
 			  
 xatom(i3)=atom(i)
 xlabel(i3)=label(i)
 xelementtype(i3)=elementtype(i)
 xamino_name(i3)=amino_name(i)
 xcol(i3)=col(i)
 xres(i3)=res(i)
 xrx(i3)=rx(i)
 xry(i3)=ry(i)
 xrz(i3)=rz(i)+boxl			  
	
	
 xatom(i4)=atom(i)
 xlabel(i4)=label(i)
 xelementtype(i4)=elementtype(i)
 xamino_name(i4)=amino_name(i)
 xcol(i4)=col(i)
 xres(i4)=res(i)
 xrx(i4)=rx(i)+boxl
 xry(i4)=ry(i)+boxl
 xrz(i4)=rz(i)	
	
 xatom(i5)=atom(i)
 xlabel(i5)=label(i)
 xelementtype(i5)=elementtype(i)
 xamino_name(i5)=amino_name(i)
 xcol(i5)=col(i)
 xres(i5)=res(i)
 xrx(i5)=rx(i)
 xry(i5)=ry(i)+boxl
 xrz(i5)=rz(i)+boxl 
 
 xatom(i6)=atom(i)
 xlabel(i6)=label(i)
 xelementtype(i6)=elementtype(i)
 xamino_name(i6)=amino_name(i)
 xcol(i6)=col(i)
 xres(i6)=res(i)
 xrx(i6)=rx(i)+boxl
 xry(i6)=ry(i)
 xrz(i6)=rz(i)+boxl
 
  xatom(i7)=atom(i)
 xlabel(i7)=label(i)
 xelementtype(i7)=elementtype(i)
 xamino_name(i7)=amino_name(i)
 xcol(i7)=col(i)
 xres(i7)=res(i)
 xrx(i7)=rx(i)+boxl
 xry(i7)=ry(i)+boxl
 xrz(i7)=rz(i)+boxl
              
			  enddo					  	  
			  
	do i=1,8*nb
		if((xrx(i).le.(boxl/2+increx)).and.(xrx(i).ge.(-boxl/2+increx)).and.(xry(i).le.(boxl/2+increy)).and.(xry(i).ge.(-boxl/2+increy)) .and.(xrz(i).le.(boxl/2+increz)).and.(xrz(i).ge.(-boxl/2+increz))) then
		   ! write(6,*)'3rd case: ',xrx(i),xry(i),xrz(i)
			check=check+1
		    if(int(MOD(i,nb)).ne.0) then 
			   m=xlabel(int(MOD(i,nb)))
			else
			!   write(6,*)'3nd case, the last bead!',j
			   m=xlabel(nb)
			endif
			
          newatom(m)=xatom(i)
		  newlabel(m)=xlabel(i)
		  newelementtype(m)=xelementtype(i)
		  newamino_name(m)=xamino_name(i)
		  newcol(m)=xcol(i)
		  newres(m)=xres(i)
		  newrx(m)=xrx(i)
	 	  newry(m)=xry(i)
		  newrz(m)=xrz(i)	
				 		
			endif
		    enddo	  		   
 !write(6,*)'check must be equal to: ',nb,check

			   endif
		        enddo
		
		close(9)

		write(6,*)'# of chain penetrating boundary: ',count
		avgnewrx=0.0
		avgnewry=0.0
		avgnewrz=0.0
		    c=0
		!center of system
		do i=1,nc
			if(count(i).eq.0) then
               c=c+1			
	           do j=1,nb
			         k=(i-1)*nb+j			
	            avgnewrx = avgnewrx+newrx(k)
                avgnewry = avgnewry+newry(k)
                avgnewrz = avgnewrz+newrz(k)	
		     enddo
		    endif
		enddo
		
		  avgnewrx = avgnewrx/c/nb
          avgnewry = avgnewry/c/nb
          avgnewrz = avgnewrz/c/nb
		!write(6,*)'physical center: ',avgnewrx,avgnewry,avgnewrz
		
!#ifdef  wholemove

!#else		
		do i=1,nc
			if(count(i).ne.0) then
			
			    do j=1,nb
			         k=(i-1)*nb+j
			avgc1x=avgc1x+newrx(k)
		    avgc1y=avgc1y+newry(k)
		    avgc1z=avgc1z+newrz(k)
				enddo
				!center of each chain
			avgc1x=avgc1x/nb
		    avgc1y=avgc1y/nb
		    avgc1z=avgc1z/nb
			
			        com1=(avgc1x-avgnewrx)**2+(avgc1y-avgnewry)**2+(avgc1z-avgnewrz)**2
					com2=(avgc1x+countx(i)*boxl-avgnewrx)**2+(avgc1y+county(i)*boxl-avgnewry)**2+(avgc1z+countz(i)*boxl-avgnewrz)**2
					com3=(avgc1x-countx(i)*boxl-avgnewrx)**2+(avgc1y-county(i)*boxl-avgnewry)**2+(avgc1z-countz(i)*boxl-avgnewrz)**2
					com4=min(com1,com2,com3)
			
					do j=1,nb
					k=(i-1)*nb+j
					
		if(com4.eq.com1) then
			newrx(k)=newrx(k)
			newry(k)=newry(k)
			newrz(k)=newrz(k)
		elseif(com4.eq.com2) then
			newrx(k)=newrx(k)+countx(i)*boxl
			newry(k)=newry(k)+county(i)*boxl
			newrz(k)=newrz(k)+countz(i)*boxl
		elseif(com4.eq.com3) then
			newrx(k)=newrx(k)-countx(i)*boxl
			newry(k)=newry(k)-county(i)*boxl
			newrz(k)=newrz(k)-countz(i)*boxl
		endif
				enddo	
				
			endif
		enddo
!#endif		

	 open(9,file=output1,action='write',form='formatted')	
		do i=1,nop
		write(9,7) newatom(i),newlabel(i),newelementtype(i),newamino_name(i),newcol(i),newres(i), newrx(i), newry(i), newrz(i)
		enddo
		
		close(9)
	 
200    write(6,*)'go to next pdb',factorx,factory,factorz
	 
       end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		end
		