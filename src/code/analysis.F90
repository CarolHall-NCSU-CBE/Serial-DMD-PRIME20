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
	INTEGER :: i,j,argcount, k,l,round, io,annealid,count,ii,jj,n, sumlist, old,iiii,m,nl,nm
	character (len = 20) :: arg(3)	
	integer izero, start, end, ioerr
	character*64 filename
	logical exist_flag
	character*1 zero
	character*1 char
	character*4 initfile, endfile,timeorcoll
	character(len=2) :: itochar
	integer ichar, len
	character*64 input1,output1,input2
	data zero/'0'/
	integer numatoms,colltmp, collcount, collanneal
	real colltotal,ttmp, annealtemp, sumanneal,presum
	real*8 ttotal,ered,tred,redtemp,realtemp,temporary,e_kin
      	real*8 rg_avg,e2e_avg
      	real*8 ehh_ii,ehh_ij
      	integer hb_alpha,hb_ii,hb_ij
	real*8 rxij,ryij,rzij,rijsq,wellsq
	integer midchain, clusternum
	integer face(6),center,sign,neighbor(50:100)
	integer, allocatable :: surface(:,:) 
	real*8 bmass_temp,bmass(28)
		
!VN: get paths for file openning and creating:
	call getcwd(rundir)
	call get_command_argument(0,path)
	realpath = scan (path,'/',back)
	mydir = path(1:realpath)
      	call readinputs()
	noptotal = nop1+nop2
	boxl = boxlength
	call allocatearrays()
	argcount = 0
  	DO
    		CALL get_command_argument(argcount, arg(argcount))
    		IF (LEN_TRIM(arg(argcount)) == 0) EXIT
		arg(argcount) = adjustl(arg(argcount))
    		!WRITE (*,*) argcount, TRIM(arg(argcount))
    		argcount = argcount+1
  	END DO

8       format(3F12.4)	

		if (arg(1) .eq. 'traj') then
			initfile = arg(2)
			endfile = arg(3)
			allocate(readposx(noptotal),readposy(noptotal),readposz(noptotal))
			call chdir(rundir)
			izero = ichar('0')
			read(initfile,*) start
			read(endfile,*) end
			round = 1
			open(999,file = 'analysis/check.txt',status = 'unknown')
			open(sf, file='analysis/run'//initfile//'to'//endfile//'.xyz',status='unknown',position='append')
			do i = start, end
				fname_digits = char(i/1000+izero)//char(mod(i,1000)/100+izero)//char(mod(i,100)/10+izero)//char(mod(i,10)+izero)
				input1 = 'results/run'//fname_digits//'.traj'
				open(traj,file = input1,status = 'old',form='unformatted')
				do while (.true.)
					
					do k = 1,noptotal
						read(traj,end=100) readposx(k), readposy(k), readposz(k)
						write(999,*) readposx(k), readposy(k), readposz(k)
					enddo
					call writesf_xyz()
				enddo
100				close(traj)
			enddo
			close(999)
!!!!!!! Energy and temperature vs time:
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
					realtemp = tred/12*2288.467-115.79
					e_kin = tred*dble(noptotal)*3.d0*0.5d0
					!!!!! Note: energies reported from PRIME20 is not in reduced unit
					write(3878,'(f15.4,f24.4,f18.1,f20.4,2f14.4)') (colltotal+real(coll)/1000000000.0d0),(ttotal+t)*0.96*0.001*3.3/sqrt(redtemp*12),realtemp,ered,e_kin,(ered-0.5*tred*3.0*noptotal)
					if((coll.gt.0).and.(mod(coll,10000000).eq.0)) goto 101
				enddo
101				close(rune)
				ttotal=ttotal+t
				colltotal = colltotal+real(coll)/1000000000.0d0
			enddo
!!!!!!! Hydrogen bonding vs time: 
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
!!!!!!! Sidechain-sidechain interaction energy vs time: 
		elseif (arg(1) .eq. 'scvst') then
			call chdir(rundir)
			initfile = arg(2)
			endfile = arg(3)
			izero = ichar('0')
			read(initfile,*) start
			read(endfile,*) end
			round = 1
			colltotal = 0.0
			ttotal = 0
			open(72, file='analysis/scvst'//initfile//'to'//endfile//'.txt',status='unknown',position='append')
			write(72,'(a25,a20,a30,a30)') 'collisions(billions)','time(microsecond)','intrapeptide sc-sc (kJ/mol)','interpeptide sc-sc (kJ/mol)'
			write(72,'(a111)') '==============================================================================================================='
			do i = start, end
				fname_digits = char(i/1000+izero)//char(mod(i,1000)/100+izero)//char(mod(i,100)/10+izero)//char(mod(i,10)+izero)
				input1 = 'results/run'//fname_digits//'.energy'
				open(rune,file = input1,status = 'old')
				if (i.gt.1) then
					read(rune,'(i15,4f12.4,3i8,4f12.4)') coll,t,redtemp,ered,tred,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij,rg_avg,e2e_avg
				endif
				do while (.true.)
					read(rune,'(i15,4f12.4,3i8,4f12.4)') coll,t,redtemp,ered,tred,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij,rg_avg,e2e_avg
					write(72,'(f15.4,2f24.4,f30.4)') (colltotal+real(coll)/1000000000.0d0),(ttotal+t)*0.96*0.001*3.3/sqrt(redtemp*12),ehh_ii,ehh_ij
					if((coll.gt.0).and.(mod(coll,10000000).eq.0)) goto 103
				enddo
103				close(rune)
				ttotal=ttotal+t
				colltotal = colltotal+real(coll)/1000000000.0d0
			enddo

!!!!!!! Intrapeptide interaction by residues: 
		elseif (arg(1) .eq. 'scvst') then
			call chdir(rundir)
			initfile = arg(2)
			endfile = arg(3)
			izero = ichar('0')
			read(initfile,*) start
			read(endfile,*) end
			round = 1
			colltotal = 0.0
			ttotal = 0
			open(72, file='analysis/scvst'//initfile//'to'//endfile//'.txt',status='unknown',position='append')
			write(72,'(a25,a20,a30,a30)') 'collisions(billions)','time(microsecond)','intrapeptide sc-sc (kJ/mol)','interpeptide sc-sc (kJ/mol)'
			write(72,'(a111)') '==============================================================================================================='
			do i = start, end
				fname_digits = char(i/1000+izero)//char(mod(i,1000)/100+izero)//char(mod(i,100)/10+izero)//char(mod(i,10)+izero)
				input1 = 'results/run'//fname_digits//'.energy'
				open(rune,file = input1,status = 'old')
				if (i.gt.1) then
					read(rune,'(i15,4f12.4,3i8,4f12.4)') coll,t,redtemp,ered,tred,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij,rg_avg,e2e_avg
				endif
				do while (.true.)
					read(rune,'(i15,4f12.4,3i8,4f12.4)') coll,t,redtemp,ered,tred,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij,rg_avg,e2e_avg
					write(72,'(f15.4,2f24.4,f30.4)') (colltotal+real(coll)/1000000000.0d0),(ttotal+t)*0.96*0.001*3.3/sqrt(redtemp*12),ehh_ii*12.47,ehh_ij*12.47
					if((coll.gt.0).and.(mod(coll,10000000).eq.0)) goto 104
				enddo
104				close(rune)
				ttotal=ttotal+t
				colltotal = colltotal+real(coll)/1000000000.0d0
			enddo

!!!!!!! Create cluster list:
		elseif (arg(1) .eq. 'cluster') then
			initfile = arg(2)
			endfile = arg(3)
			izero = ichar('0')
			read(initfile,*) start
			read(endfile,*) end
			call chdir(rundir)
			
			open(258, file='analysis/clustervst'//initfile//'to'//endfile//'.txt',status='unknown',position='append')
			write(258,'(a25,a20,a11,a11,a9)') 'collisions(billions)','time(microsecond)','Monomer','Oligomer','Fibril'
			write(258,'(a88)') '================================================================================================'
			do i = start, end
				call chdir(rundir)
				fname_digits = char(i/1000+izero)//char(mod(i,1000)/100+izero)//char(mod(i,100)/10+izero)//char(mod(i,10)+izero)
				input1 = 'results/run'//fname_digits//'.energy'
				
				open(rune,file = input1,status = 'old')
				if (i.gt.1) then
					read(rune,'(i15,4f12.4,3i8,4f12.4)') coll,t,redtemp,ered,tred,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij,rg_avg,e2e_avg
				endif
				do while (.true.)
					read(rune,'(i15,4f12.4,3i8,4f12.4)') coll,t,redtemp,ered,tred,hb_alpha,hb_ii,hb_ij,ehh_ii,ehh_ij,rg_avg,e2e_avg
					collinbill = colltotal+real(coll)/1000000000.0d0
					write(collnum,'(f7.4)') collinbill
					call readconfig
					call readbptnr
					call readparameters
					call clustertrack
					write(258,'(f15.4,f24.4,i14,2i10)') (colltotal+real(coll)/1000000000.0d0),(ttotal+t)*0.96*0.001*3.3/sqrt(redtemp*12),coil,oligomer,fibril
					if((coll.gt.0).and.(mod(coll,10000000).eq.0)) goto 107
				enddo
107				close(rune)
				ttotal=ttotal+t
				colltotal = colltotal+real(coll)/1000000000.0d0
			enddo
			do i = 1, (nc+nc2)
				write(258,'(50i5)') pack(layer(i,:), layer(i,:)/= 0)
			enddo
					do i = 1,(nc+nc2)
						res  = 0
						call phipsianalysis(i)
						do l=0,360
         						do m=0,360
	    							if (res(l,m) .ne. 0) then
               								nl=l-180
               								nm=m-180
									!write(258,*) i,nl,nm
	    							endif 
        						enddo
						enddo
					enddo
								
			do i = 1, (noptotal-1)
				do j = (i+1), noptotal
					if ((bptnr(i) .eq. j) .and. (chnnum(i) .ne.chnnum(j))) then
						!write(258,*) chnnum(i), chnnum(j), i,j
					endif
				enddo
			enddo	
		close(258)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		elseif (arg(1) .eq. 'aggregation') then
			call chdir(rundir)
			initfile = arg(2)
			endfile = arg(3)
			izero = ichar('0')
			read(initfile,*) start
			read(endfile,*) end
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
					write(3878,'(f15.4,f24.4,f18.1,f20.4,2f14.4)') (colltotal+real(coll)/1000000000.0d0),(ttotal+t)*0.96*0.001*3.3/sqrt(redtemp*12),realtemp,ered*12.47,tred*12.47,(ered-0.5*tred*3.0*noptotal)*12.47
					if((coll.gt.0).and.(mod(coll,10000000).eq.0)) goto 501
				enddo
501				close(rune)
				ttotal=ttotal+t
				colltotal = colltotal+real(coll)/1000000000.0d0
			enddo



		
!!!!!!! Write pdb file at any collision:
		elseif (arg(1) .eq. 'pdb') then
			call chdir(rundir)
			collnum = arg(2)
			read(collnum,*) collinbill
			call readparameters
			call readconfig
			call readbptnr
			! Remove periodic boundary condition:
			call clustertrack
			open(runpdb, file = 'analysis/'//collnum//'billioncollision.pdb',status = 'unknown')
			do k=1,(noptotal)
				sv(1,k) = sv(1,k)*boxl
				sv(2,k) = sv(2,k)*boxl
				sv(3,k) = sv(3,k)*boxl
      			enddo				

			do k = 1, nop1
	 			chnnum(k)=(k-1)/numbeads1+1
      			end do
      			do k = 1, nop2
	 			chnnum(nop1+k)=(nop1/numbeads1)+(k-1)/numbeads2+1
      			end do

			call removepbc
			call write_rasmol

!!!!!!! Test bptnr:
		elseif (arg(1) .eq. 'test') then
		call chdir(rundir)
		open(333,file='analysis/test.txt',status = 'unknown',position = 'append')
		open(unit=111,file='results/run0010.config',status = 'old',form='unformatted')
		open(unit=11,file='results/run0010.bptnr',status = 'old',form='unformatted')
		do while (.true.)
			read(111) coll,t,old_rx,old_ry,old_rz
			if((coll.gt.0).and.(mod(coll,10000000).eq.0)) goto 222
		enddo
222	close(111)	
		
		do while (.true.)
			read(11) coll,bptnr
			if((coll.gt.0).and.(mod(coll,10000000).eq.0)) goto 223
		enddo
223	close(11)
		do k = 1, nop1
	 			chnnum(k)=(k-1)/numbeads1+1
      			end do
      			do k = 1, nop2
	 			chnnum(nop1+k)=(nop1/numbeads1)+(k-1)/numbeads2+1
      			end do

		do i = 1,noptotal
			sv(1,i) = old_rx(i)
			sv(2,i) = old_ry(i)
			sv(3,i) = old_rz(i)
			write(333,*) i,bptnr(i),chnnum(i),chnnum(bptnr(i))
		enddo
		
		do j = 1,(nc+nc2)
			res  = 0
			call phipsianalysis(j)
			do l=0,360
         			do m=0,360
	    				if (res(l,m) .ne. 0) then
   						nl=l-180
               					nm=m-180
						write(333,*) j,nl,nm
	    				endif 
        			enddo
			enddo
		enddo

	endif	
	end program
#include "readinputs.f"
#include "writesf_xyz.f"
#include "allocatearrays.f"
#include "write_rasmol-YM.f"
#include "readconfig.f"
#include "readparameters.f"
#include "readbptnr.f"
#include "cluster.f"
#include "phipsi_analysis.f"
#include "removepbc.f"
