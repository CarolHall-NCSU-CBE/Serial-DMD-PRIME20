! =========================================================================
! This program read inputs from text file and generate initial configuration
! This is a serial job
! Last modified on 10/12/2023 by Van Nguyen
! =========================================================================  
#include "def.h"
#include "header.f"
 
      	program initconfig
        
      	use global
	use inputreadin

      	implicit none

      	logical over,update,xpulse_del
      	character*64 filename
	real, allocatable :: annealtemp(:)  
      	real*8 tarray_nv(2),delta_nv,extime_nv,colrat_nv
      	real*8 tgho,xr,rating,avegtime,ered,tred,sumvel,e_int
      	real*8 old_tfalse,w,v1,v2,v3,r,fact,ran_non,ran_brk,tstar
      	real*8 rxij,ryij,rzij,vxij,vyij,vzij,bij,vijsq,rijsq,diff,rg_avg,e2e_avg
	integer totalsteps
	real*8 start, finish
	character(len=1) :: c

#ifdef equil
      	real per_hb,per_hh
#endif

#ifndef canon
      	real*8 pre_energy,e_pot,e_kin
#endif

      	integer boundbad, unboundbad, pre_boundbad,pre_unboundbad,hb_partner
      	integer i,j,ii,jj,k,k_j,kk,kk_j,kkk,l,m,mm,num_before_sheet,evcode,writeprop_count
      	integer nbrsum,nl,nm,hb_alpha,hb_ii,hb_ij,hh_ii,hh_ij,n_interval
      	integer*8 nupdates,numghosts,nforcedupdate
      	integer*8 ncoll_max,nevents(30)
      	real*8  ehh_ii,ehh_ij

#ifdef write_phipsi
      	integer nphipsi
#endif

#ifdef __ifc
      integer ifc_time(2),idumb
#endif

!VN: get paths for file openning and creating:
	call getcwd(rundir)
	call get_command_argument(0,path)
	realpath = scan (path,'/',back)
	mydir = path(1:realpath)
      	data nevents / 30*0/
      	pi=4.d0*datan(1.d0)
	call readinputs()
	noptotal = nop1+nop2
	open(unit=1717,file='inputs/simtemp',status='new')
	write(1717,*) simtemp
	write(1717,*) simcoll
	close(1717)

	if (annealcheck == 0) then
		annealingsteps = 9
		allocate(annealtemp(annealingsteps))
		annealtemp= (/0.50,0.45,0.40,0.35,0.30,0.28,0.26,0.24,0.22/)	
		allocate(collset(annealingsteps))
		collset(1:5) = 100000000
		collset(6:7) = 150000000
		collset(8) = 200000000
		collset(9) = 250000000
		do j = 1,9
			write(c, '(i0)') j
			open(unit=1717,file='inputs/annealtemp_'//c,status='new')
			write(1717,*) annealtemp(j)
			write(1717,'(i9)') collset(j)
			close(1717)
		enddo
	elseif (annealcheck == 1) then
		annealingsteps = (maxannealing-minannealing)/incrannealing+1
		allocate(annealtemp(annealingsteps))
		do j = 1, annealingsteps
			annealtemp(j) = ((maxannealing -(j-1)*incrannealing)+115.79)/2288.46
		enddo
		allocate(collset(1))
		collset = annealcoll
		do j = 1,annealingsteps
			write(c, '(i0)') j
			open(unit=1717,file='inputs/annealtemp_'//c,status='new')
			write(1717,*) annealtemp(j)
			write(1717,'(i9)') collset
			close(1717)
		enddo
	endif
	call allocatearrays()
	call genconfig()
						
	End	
#include "readinputs.f"
#include "genconfig.f"
#include "bumped.f"      
#include "sqshlder.f"      
#include "check_nc_int.f"
#include "add_tbin.f"
#include "bond.f"
#include "cell_add.f"
#include "cell_link.f"
#include "check_sigma.f" 
#include "checkover.f"
#include "config.f"
#include "core.f"
#include "displ.f"
#include "energy.f"
#include "eventdyn.f"
#include "eventredo_up.f"
#include "eventredo_down.f"
#include "events.f"
#include "del_tbin.f"
#include "inputinfo.f"
#include "make_code.f"
#include "nbor.f"
#include "nbor_setup.f"
#include "nc_sqwel.f"
#include "partial_events.f"
#include "files_opn.f"
#include "write_rasmol-YM.f"
#include "repuls_add.f"
#include "repuls_check.f"
#include "repuls_check_3.f"
#include "repuls_del_a.f"
#include "repuls_del_b.f"
#include "scale_down.f"
#include "scale_up.f"
#include "sqwel.f"
#include "files_close.f"
#include "allocatearrays.f"
!LR: I separated the analysis code out of the main code.
#include "radgyr.f"  
#include "end2end.f"

#ifdef equil
#include "equil.f"
#endif

#ifdef write_phipsi
#include "phipsi.f"
#endif
