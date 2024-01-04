! =========================================================================
! This program read inputs from text file and analyse results
! This is a serial job
! Last modified on 11/14/2023 by Van Nguyen
! =========================================================================  
#include "def.h"
#include "header.f"
 
      	program dmdanalysis
        
      	use global
	use inputreadin

      	implicit none
	INTEGER :: i,j,argcount, k,round, io
  	!CHARACTER(len=32) :: arg()
	character (len = 5) :: arg(3)	
	integer izero, start, end
	character*64 filename
	logical exist_flag
	character*1 zero
	character*1 char
	character*4 initfile, endfile
	integer ichar, len
	character*64 input1,output1,input2
	!real*4, allocatable :: readpos(:,:)
	data zero/'0'/
	integer numatoms
	real colltotal
	real*8 ttotal,ered,tred,redtemp,realtemp
      	real*8 rg_avg,e2e_avg
      	real*8 ehh_ii,ehh_ij
      	integer hb_alpha,hb_ii,hb_ij
	

!VN: get paths for file openning and creating:
	call getcwd(rundir)
	call get_command_argument(0,path)
	realpath = scan (path,'/',back)
	mydir = path(1:realpath)
      	call readinputs()
	call allocatearrays()
	noptotal = nop1+nop2
	allocate(readposx(noptotal),readposy(noptotal),readposz(noptotal))
	argcount = 0
  	DO
    		CALL get_command_argument(argcount, arg(argcount))
    		IF (LEN_TRIM(arg(argcount)) == 0) EXIT
		arg(argcount) = adjustl(arg(argcount))
    		!WRITE (*,*) argcount, TRIM(arg(argcount))
    		argcount = argcount+1
  	END DO

8       format(3F12.4)	
	!do j = 1, (argcount-1)
		if (arg(1) .eq. 'traj') then
			initfile = arg(2)
			endfile = arg(3)
			call chdir(rundir)
			izero = ichar('0')
			read(initfile,*) start
			read(endfile,*) end
			round = 1
			open(sf, file='analysis/run'//initfile//'to'//endfile//'.xyz',status='unknown',position='append')
			do i = start, end
				fname_digits = char(i/1000+izero)//char(mod(i,1000)/100+izero)//char(mod(i,100)/10+izero)//char(mod(i,10)+izero)
				input1 = 'results/run'//fname_digits//'.traj'
				open(traj,file = input1,status = 'old',form='unformatted')
				!write(*,*) input1
				do while (.true.)
					do k = 1,noptotal
					read(traj,end=100) readposx(k), readposy(k), readposz(k)
					!write(999,*) readposx(k), readposy(k), readposz(k)
					enddo
					call writesf_xyz()
				enddo
100				close(traj)
			enddo
			close(999)
		elseif (arg(1) .eq. 'evst') then
			call chdir(rundir)
			initfile = arg(2)
			endfile = arg(3)
			izero = ichar('0')
			read(initfile,*) start
			read(endfile,*) end
			round = 1
			colltotal = 0.0
			ttotal = 0
			open(3878, file='analysis/evst'//initfile//'to'//endfile//'.txt',status='unknown',position='append')
			write(3878,'(a25,a20,a17,a17,a13,a13)') 'collisions(billions)','time(microsecond)','temperature(K)','Etotal(kJ/mol)','KE(kJ/mol)','PE(kJ/mol)'
			write(3878,'(a111)') '================================================================================================================'
			do i = start, end
				fname_digits = char(i/1000+izero)//char(mod(i,1000)/100+izero)//char(mod(i,100)/10+izero)//char(mod(i,10)+izero)
				input1 = 'results/run'//fname_digits//'.energy'
				open(rune,file = input1,status = 'old')
				if (i.gt.1) then
					read(rune,'(i15,4f12.4,3i8,4f12.4)') coll,t,redtemp,ered,tred,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij,rg_avg,e2e_avg
				endif
				do while (.true.)
					read(rune,'(i15,4f12.4,3i8,4f12.4)') coll,t,redtemp,ered,tred,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij,rg_avg,e2e_avg
					realtemp = redtemp*2288.467-115.79
					write(3878,'(f15.4,f24.4,f18.1,f20.2,2f14.2)') (colltotal+real(coll)/1000000000.0d0),(ttotal+t)*0.96*0.001*3.3/sqrt(redtemp*12),realtemp,ered*12.47,tred*12.47,(ered-0.5*tred*3*noptotal)
					if((coll.gt.0).and.(mod(coll,10000000).eq.0)) goto 101
				enddo
101				close(rune)
				ttotal=ttotal+t
				colltotal = colltotal+real(coll)/1000000000.0d0
			enddo
		elseif (arg(1) .eq. 'hbvst') then
			call chdir(rundir)
			initfile = arg(2)
			endfile = arg(3)
			izero = ichar('0')
			read(initfile,*) start
			read(endfile,*) end
			round = 1
			colltotal = 0.0
			ttotal = 0
			open(428, file='analysis/hbvst'//initfile//'to'//endfile//'.txt',status='unknown',position='append')
			write(428,'(a25,a20,a11,a19)') 'collisions(billions)','time(microsecond)','total hb','interpeptide hb'
			write(428,'(a77)') '============================================================================'
			do i = start, end
				fname_digits = char(i/1000+izero)//char(mod(i,1000)/100+izero)//char(mod(i,100)/10+izero)//char(mod(i,10)+izero)
				input1 = 'results/run'//fname_digits//'.energy'
				open(rune,file = input1,status = 'old')
				if (i.gt.1) then
					read(rune,'(i15,4f12.4,3i8,4f12.4)') coll,t,redtemp,ered,tred,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij,rg_avg,e2e_avg
				endif
				do while (.true.)
					read(rune,'(i15,4f12.4,3i8,4f12.4)') coll,t,redtemp,ered,tred,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij,rg_avg,e2e_avg
					write(428,'(f15.4,f24.4,2i15)') (colltotal+real(coll)/1000000000.0d0),(ttotal+t)*0.96*0.001*3.3/sqrt(redtemp*12),(hb_ii+hb_ij),hb_ij
					if((coll.gt.0).and.(mod(coll,10000000).eq.0)) goto 102
				enddo
102				close(rune)
				ttotal=ttotal+t
				colltotal = colltotal+real(coll)/1000000000.0d0
			enddo
		endif
	!enddo
		
	end program
#include "readinputs.f"
#include "writesf_xyz.f"
#include "allocatearrays.f"
