        subroutine hb_dssp(hb_alpha)

#include "def.h"

        USE GLOBAL
	use inputreadin
                      
	IMPLICIT NONE

        !integer i,j,chainn,n,ii,EVCODE
	real*8 RXIJ,RYIJ,RZIJ,RIJSQ,BLMAXSQ,rx1,ry1,rz1,r1
        character*1 col(nop1/numbeads1+nop2/numbeads2)
        character*3 amino_name
        integer i,j,k,m,n,p
	integer atboundary(3,noptotal), clusteratpbc(nc+nc2),sumatpbc(6)
	integer atpbc(2), chnatpbc((nc+nc2),6), count,clustersize,pos,mainface,nummove,sumchain
	integer hb_alpha
	call chdir(rundir)
	!!!!!!! alpha_helix (intrapeptide)
	hb_alpha = 0
	do k=1,(nop1/numbeads1)
		do i=(k-1)*numbeads1+chnln1+5,(k-1)*numbeads1+2*chnln1
			j=i+chnln1-4
			if (bptnr(i) == j) hb_alpha = hb_alpha + 1
		enddo
	enddo

	do k=1,(nop2/numbeads2)
		do i=nop1+(k-1)*numbeads2+chnln2+5,nop1+(k-1)*numbeads2+2*chnln2
			j=i+chnln2-4
			if (bptnr(i) == j) hb_alpha = hb_alpha + 1
		enddo
	enddo

	!!!!!!! Bridge (interpeptide)
	return
	end		