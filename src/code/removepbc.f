        subroutine removepbc

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
	integer face(6)
	call chdir(rundir)
	!open(8367,file='analysis/temp.txt',status='unknown')
	atboundary = 0
	chnatpbc = 0
	clusteratpbc = 0
	count = 0
	face = 0
	atpbc = maxloc(abs(sv(1:3,:)))	
	!!!!!! Finding main face = node 0
	do i = 1,noptotal
		do j = 1,3
			if ((boxl/2-abs(sv(j,i))) .lt.1.0) then
				if (sv(j,i) .lt. 0) then
					chnatpbc(chnnum(i),j) = 1
					!print*, j, chnnum(i)
				else
					chnatpbc(chnnum(i),j+3) = 1
					!print*, j+3, chnnum(i)
				endif
			endif
		enddo
	enddo
	sumatpbc = sum(chnatpbc,dim=1)
	!print*, sumatpbc
	nummove = size(pack(sumatpbc, sumatpbc/= 0))/2
	!print*, nummove
	!!!!!!! Identify how many clusters being segmented by pbc
	do i = 1,(nc+nc2)
		sumchain = 0
		do j = 1, (nc+nc2)
			!write(8367,*) i, j, cluster(i,j)
			if (cluster(i,j).ne.0) then
				sumchain = sumchain + sum(chnatpbc(cluster(i,j),:))
			endif
		enddo
		if (sumchain.ne.0) then	
			count = count + 1
			!print*, i, count
			clusteratpbc(count) = i
		endif
		!write(8367,'(6i6)') chnatpbc(i,:)
		
	enddo
	do i = 1,count
		pos = clusteratpbc(i)
		!print*, 'count: ', i, 'pos: ', pos
		clustersize = size(pack(cluster(pos,:), cluster(pos,:)/= 0))
		!print*, 'clustersize: ', clustersize
		do j = 1,(nc+nc2)
			if ((cluster(pos,j).ne.0).and.(cluster(pos,j) .le. nc)) then
				!print*, 'chain in cluster: ', cluster(pos,j)
				do k = chnln1*(cluster(pos,j)-1)+1, chnln1*cluster(pos,j)
					do m = 1,3
						if((sv(m,k).gt.(boxl/10*9)).and.(sumatpbc(m).ne.0)) then
							face(m) = face(m)+1
						elseif ((sv(m,k).lt.(-boxl/10*9)).and.(sumatpbc(m+3).ne.0)) then
							face(m+3) = face(m+3)+1
						endif
					enddo
				enddo
			elseif ((cluster(pos,j).ne.0).and.(cluster(pos,j) .gt. nc)) then
				!print*, 'chain in cluster: ', cluster(pos,j)
				do k =  chnln2*(cluster(pos,j)-1)+1, chnln2*cluster(pos,j)
					do m = 1,3
						if((sv(m,k).gt.0).and.(sumatpbc(m).ne.0)) then
							face(m) = face(m)+1
						elseif ((sv(m,k).lt.0).and.(sumatpbc(m+3).ne.0)) then
							face(m+3) = face(m+3)+1
						endif
					enddo
				enddo
			endif
		enddo
		!!!!!!! Finding edges between nodes - shortest distance to move cluster
		do n = 1,nummove
			mainface = maxloc(face,dim=1)
			!print*, mainface
			do j = 1,(nc+nc2)
				if ((cluster(pos,j).ne.0).and.(cluster(pos,j) .le. nc)) then
					do k = numbeads1*(cluster(pos,j)-1)+1, numbeads1*cluster(pos,j)
						if (((int(real(mainface)/4)*(-2)+1) *sv((int(real(mainface)/4)*(-3))+mainface,k)).lt.0) then
							sv((int(real(mainface)/4)*(-3))+mainface,k) = sv((int(real(mainface)/4)*(-3))+mainface,k) + (int(real(mainface)/4)*(-2)+1)*boxl
							
					endif
					enddo
				elseif ((cluster(pos,j).ne.0).and.(cluster(pos,j) .gt. nc)) then
					do k =  numbeads2*(cluster(pos,j)-1)+1, numbeads2*cluster(pos,j)
						if (((int(real(mainface)/4)*(-2)+1)*sv((int(real(mainface)/4)*(-3))+mainface,k)).lt.0) then
							sv((int(real(mainface)/4)*(-3))+mainface,k) = sv((int(real(mainface)/4)*(-3))+mainface,k) + (int(real(mainface)/4)*(-2)+1)*boxl
						endif
					enddo
				endif
			enddo
			face(mainface) = 0
		enddo			
	enddo		
							
	!close(8367)
	
	end		