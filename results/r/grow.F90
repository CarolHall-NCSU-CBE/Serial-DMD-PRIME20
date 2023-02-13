#include "def.h"
#include "header.f"

! % of peptide decrease or increase

        program analysis_code
        use global
        implicit none  
		
          integer i,j,k,kk,iii,ii,jj,l,m,n,count_file2
		  integer count,recordid(13)
          logical start_file
        character*64 filename
        character*4 fname_digits
        common fname_digits		
        integer izero,iflag
        logical exist_flag
        character*1 zero,char
        integer ichar
        data zero/'0'/
	
        izero = ichar('0')
        iflag = 9
        filename = '../run'//'0009.config'
        fname_digits = '0009'
        inquire(file = filename, exist = exist_flag)
         ttotal = 0.0
         t = 0.0
8       format(A4,3X,I4,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3)		
2 	    format(i15,3f12.4,3i8,4f12.4)	

        start_file = .false.
         pi=4.d0*datan(1.d0)
		 ncoll=1000000000
		noptotal=nop1+nop2
		count_file = 0
        hb_total = 0
		std1 = 0.0 
		hb2_total = 0 
		std2 = 0.0 
		hb_cross_total = 0
		std3 = 0.0
       count = 0
		!boxl=110.0d0

      	open(unit=7,file='../../parameters/identity.inp',status='unknown')
      	read (7,*) ab
      	close(unit=7)
		
      	do l=1,nop1,numbeads1
		    do k=1,numbeads1
			  identity(l+k-1)=ab(k)
		    enddo
	      enddo		
      	do l=nop1+1,nop1+nop2,numbeads2
         	do k=1,numbeads2
               identity(l+k-1)=ab(numbeads1+k)
         	enddo
      	enddo	  
	  
      	open(unit=7,file='../../parameters/firstside1.data',status='unknown')
      	read(7,*) fside1
      	close(unit=7)
      	l=1
      	do k=1,chnln1
         	if (fside1(k) .ne. 0) then
	    		fside1(k)=3*chnln1+l
	    		l=l+1
	 	    endif
      	enddo	
! from no Gly to with Gly
      	open(unit=7,file='../../parameters/firstside2.data',status='unknown')
      	read(7,*) fside2
      	close(unit=7)
      	l=1
      	do k=1,chnln2
            fside2(k)=3*chnln2+l
	      l=l+1	
      	enddo	
		
      	do k = 1,nop1
	 	    chnnum(k)=(k-1)/numbeads1+1
      	end do
      	do k = 1,nop2
	 	    chnnum(nop1+k)=nc1+(k-1)/numbeads2+1	
      	end do

      	do k = 1,nc1
	      do kk=1,chnln1     !from 1-40
		    i = (k-1)*numbeads1+kk
	 	    residue(i)=kk
		    ii = (k-1)*numbeads1+chnln1+kk
	 	    residue(ii)=kk
		    iii = (k-1)*numbeads1+chnln1*2+kk
	 	    residue(iii)=kk
	      enddo
      	enddo		
		
      	do k = 1,nc2
	      do kk=1,chnln2  ! from 1-7
		    i = nop1+(k-1)*numbeads2+kk
	 	    residue(i)=kk
		    ii = nop1+(k-1)*numbeads2+chnln2+kk
	 	    residue(ii)=kk
		    iii = nop1+(k-1)*numbeads2+chnln2*2+kk
	 	    residue(iii)=kk
	      enddo
      	enddo			

      	open(unit=7,file='../../parameters/beadwell_ha55a.data',status='unknown')
      	do i=1,(20)*(20)
      		read(7,711) ii,jj,bdtemp,wltemp
711  		format(2(2x,i2,2x),2(f6.3,2x))
       		bds(ii,jj)=bdtemp
       		wel(ii,jj)=wltemp
      	enddo
	     close(unit=7)
!       open(unit=7,file='../../parameters/beadwell_ha45a.data',status='unknown')
!        do i=1,400
!           read(7,711) ii,jj,bdtemp,wltemp
!           wel2(ii,jj)=wltemp
!        end do
!        close(unit=7)

      	open(7,file='../../parametersep/ep19p_ha55a_weakhp.data',status='unknown')
      	do i=1,(20)*(20)
      		read(7,"(2(2x,i2,2x),f8.3)") ii,jj,eptemp
      		ep(ii,jj)=eptemp!*0.7
		!	ep2(ii,jj)=eptemp*1.3
      	end do
      	close(7)		
		
      open(9,file='data.dat',action='write',form='formatted')
         open(10,file='oligomer.dat',action='write',form='formatted')  
	  
       do while (exist_flag)
         !  if (iflag.ge.stopfile) goto 333 	   
		iflag = iflag + 1
		fname_digits = char(iflag/1000+izero)//char(mod(iflag,1000)/100+izero)//char(mod(iflag,100)/10+izero)//char(mod(iflag,10)+izero)
        inquire( file = '../run'//fname_digits//'.pdb', exist = exist_flag)

      if (exist_flag) then
             else 
          write(6,*)'reach end file minus 2',iflag
           goto 333
       endif	 
		  count_file2  = 0
	     open(7,file='../run'//fname_digits//'.config', form='unformatted')
         open(8,file='../run'//fname_digits//'.bptnr', form='unformatted')  
		 
		 ttotal = ttotal + t_temp
          t_temp = 0.0
          do while (.true.)			 
     	     read(8,end=120) coll,bptnr
       		 read(7,end=120) coll,t,rx,ry,rz
    count_file2 = count_file2 + 1		 
			 
       ! if ((count_file2.eq.1).and.(fname_digits.eq.'0175')) then
   !if (iflag.ge.610) then
           call fibril_list_assign()  ! Yiming

      write(9,"(f12.3,2x,a4,29(2x,i4))") ttotal+t,fname_digits,target_pep_num_in_fibril,sheet_num,fibril_num,num_sheet(target_fibril),&
&  pep_num(1),pep_num(2),pep_num(3),pep_num(4),pep_num(5),pep_num(6),pep_num(7),pep_num(8),pep_num(9),pep_num(10),pep_num(11),pep_num(12),pep_num(13),pep_num(14),&
&  pep_num(15),pep_num(16),pep_num(17),pep_num(18),pep_num(19),pep_num(20),pep_num(21),pep_num(22),pep_num(23),pep_num(24),pep_num(25)
         !  write(6,*) 'peptides within fibril ',ttotal+t,target_pep_num_in_fibril
         !  start_file = .true.
    ! endif	
    	  if (fname_digits.eq.'0093') then
            do i=1,sheet_num
            do j= 1,nc1+nc2
			!	 write(6,*) j,sheet_identity(j)
			  if(sheet_identity(j).eq.i) then
                count = count + 1			
			    recordid(count) = j
			  endif
		  enddo	
		        write(10,"(15(2x,i4))") i,pep_num(i),recordid(1),recordid(2),recordid(3),recordid(4),recordid(5),recordid(6),recordid(7),recordid(8),recordid(9),recordid(10),recordid(11),recordid(12),recordid(13)
	       
		   count = 0
		   do j=1,13
		   recordid(j) = 0
		   enddo
		   
 		   enddo
        endif	
			   
	        count_file = count_file + 1  
		   
              if(t_temp.le.t) t_temp = t 
 
         enddo	  	
120   continue	       

	        close(7)
		    close(8)
  ! reading next file
        end do
					 
333	      close(9)		
          close(10)			  
		    end
		
#include "fibril_list_assign.f"
!#include "fibril_list_update.f"