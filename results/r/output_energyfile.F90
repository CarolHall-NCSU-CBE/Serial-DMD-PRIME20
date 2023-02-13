        program fibrilkinetics
		
       implicit none

	
       integer, parameter :: k12=selected_int_kind(12)		
       integer(kind=k12) coll,ncoll		
		real t,ttotal
		integer noptotal
      	real*8 ered,tred
      	real*8 rg_avg,e2e_avg
      	real*8  ehh_ii,ehh_ij
      	integer hb_alpha,hb_ii,hb_ij,countattach(500),countfile
        integer check,a,b,c,i,j,k,l,m,n,p,q,o

!		real*8 rx(nop),ry(nop),rz(nop)
 !     integer, dimension(:), allocatable::cell,wrap_map			 

		logical attach,skip

		character*64 input1,output1,input2
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
	
 !     	allocate(fibrilnop(nc1*nb))
		
8       format(A4,3X,I4,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3)		
2 	    format(i15,3f12.4,3i8,4f12.4)	
		ncoll=1000000000
		noptotal=64
			!boxl=110.0d0
			
	 open(9,file='polyala_hb.dat',action='write',form='formatted')				
	! open(11,file='.dat',action='write',form='formatted')		 

	do while (exist_flag)		
		iflag = iflag + 1

		fname_digits = char(iflag/1000+izero)//char(mod(iflag,1000)/100+izero)//char(mod(iflag,100)/10+izero)//char(mod(iflag,10)+izero)
	    input1='../run'//fname_digits//'.pdb'

	   ! if(iflag.eq.501) goto 777
	 if (exist_flag) then
	 else 
	 write(6,*)'reach end file minus 2',iflag
	 goto 333
    endif	

	    open(7,file='../run'//fname_digits//'.energy',status='unknown')
      	do while (.true.)
      	read(7,2) coll,t,ered,tred,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij,rg_avg,e2e_avg
		write(9,*) ttotal+t, hb_alpha,hb_ii-hb_alpha,hb_ij
		  if((coll.gt.0).and.(mod(coll,10000000).eq.0)) goto 111
       	enddo

111		   close(7)		
		
!      	deallocate(fibrilnop)	 
         
	      ttotal=ttotal+t
          
			 
  ! reading next file

       end do
333	      close(9)
		  !  close(11)		
		    end
		
