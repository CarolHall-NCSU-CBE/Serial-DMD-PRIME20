
   
#include "def.h"
#include "header.f"
 
      	program dmd

      	use global
	use inputreadin
	

      	implicit none

      	logical over,update,xpulse_del
      	character*64 filename  
      	real*8 tarray_nv(2),delta_nv,extime_nv,colrat_nv
      	real*8 tgho,xr,rating,avegtime,ered,tred,sumvel,e_int
      	real*8 old_tfalse,w,v1,v2,v3,r,fact,ran_non,ran_brk,tstar
      	real*8 rxij,ryij,rzij,vxij,vyij,vzij,bij,vijsq,rijsq,diff,rg_avg,e2e_avg

#ifdef equil
      	real per_hb,per_hh
#endif

#ifndef canon
      	real*8 pre_energy,e_pot,e_kin
#endif

      	integer boundbad, unboundbad, pre_boundbad,pre_unboundbad,hb_partner
      	integer i,j,ii,jj,k,k_j,kk,kk_j,l,m,mm,num_before_sheet,evcode,writeprop_count
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

!LR: Calculates the total number of particles in the system.  I changed this to be a modular variable,
! as there are a number of places that this streamlines the code, and it will hopefully prevent some of the 
! tediousness that I had in adapting Dave's non-modular 2-species code to a modular 3-species code.
	
	!VN: get paths for files:
	call getcwd(rundir)
	call get_command_argument(0,path)
	realpath = scan (path,'/',back)
	mydir = path(1:realpath)

end