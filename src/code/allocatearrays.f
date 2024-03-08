	subroutine allocatearrays()
   
#include "def.h"

	use global
	use inputreadin
 
	implicit none
	
	allocate (old_rx(noptotal),old_ry(noptotal),old_rz(noptotal))
	allocate (sv(6,noptotal))
	allocate (ev_code(noptotal,noptotal))
	allocate (hp(numbeads1+numbeads2),hp1(numbeads1),hp2(numbeads2))
	allocate (chnnum(noptotal),coltype(noptotal+3),extra_repuls(noptotal,4))
	allocate (bptnr(noptotal),identity(noptotal),nptnr(noptotal+3))
	allocate (npt(noptotal+1),nb(maxnbs*noptotal+1),na_npt(noptotal))
	allocate (npt_dn(noptotal+1),dnnab(maxnbs*noptotal+1),nnabdn(noptotal))
	allocate (tlinks(noptotal+3),tlinks2(noptotal+3))
	allocate (clinks(noptotal),fside1(chnln1),fside2(chnln2))
	allocate (bdln(chnln1+chnln2),bl_rn(chnln1+chnln2),bl_rc(chnln1+chnln2))
	allocate (tim(noptotal+3), bm(noptotal))
	allocate (del_bdln(chnln1+chnln2),del_blrn(chnln1+chnln2),del_blrc(chnln1+chnln2))
	allocate (cluster(nc+nc2,nc+nc2),layer(nc+nc2,nc+nc2))
	
      	return

      	end




