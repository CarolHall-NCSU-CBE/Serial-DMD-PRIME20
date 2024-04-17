subroutine clustertrack
	use global
	use inputreadin
	implicit none
	
	!!!!!!! Variable for secondary structure identification:
	integer connect
	integer, allocatable :: hblist(:,:),hplist(:,:),link(:),linklayer(:)
	INTEGER :: i,j,k,l,ii,jj,count,clustersize,m,nl,nm,bres,layersize
	real*8 rxij,ryij,rzij,rijsq,wellsq


	
	allocate(hblist(nc+nc2,nc+nc2),hplist(nc+nc2,nc+nc2),link(nc+nc2),linklayer(nc+nc2))
	call chdir(rundir)
			hblist = 0
			hplist = 0
			cluster = 0
			layer = 0
			link = 0
			coil = 0
			dimer = 0
			oligomer = 0
			fibril = 0
			do k = 1, nop1
	 			chnnum(k)=(k-1)/numbeads1+1
      			end do
      			do k = 1, nop2
	 			chnnum(nop1+k)=(nop1/numbeads1)+(k-1)/numbeads2+1
      			end do
			
			!!!!!!! Counting Hydrogen bonds between peptides
			do i = 1, (noptotal-1)
				do j = (i+1), noptotal
					if ((bptnr(i) .eq. j) .and. (chnnum(i) .ne.chnnum(j))) then
						hblist(chnnum(i),chnnum(j)) = hblist(chnnum(i),chnnum(j)) + 1
					endif
					if ((identity(i).gt.8).and.(identity(j).gt.8)) then
						if (i .le. nop1) then
							ii=i-(chnnum(i)-1)*numbeads1
						elseif (i .le. nop1+nop2) then
          						ii=i-nop1-((chnnum(i)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
          					endif
	    				if (j .le. nop1) then
	    					jj=j-(chnnum(j)-1)*numbeads1
           				elseif (j .le. nop1+nop2) then
           					jj=j-nop1-((chnnum(j)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
           				endif
                 				if ((hp(ii).eq.1).and.(hp(jj).eq.1)) then
                    					if ((chnnum(i).ne.chnnum(j))) then
                       						rxij=sv(1,i)-sv(1,j)
               							ryij=sv(2,i)-sv(2,j)
               							rzij=sv(3,i)-sv(3,j)
               							rxij=rxij-dnint(rxij)
               							ryij=ryij-dnint(ryij)
               							rzij=rzij-dnint(rzij)
               							rijsq=rxij*rxij+ryij*ryij+rzij*rzij
								wellsq=welldia_sq(identity(i),identity(j))
								!print*, i, j,identity(i),identity(j), chnnum(i), chnnum(j), wellsq
               							if (rijsq.le.wellsq) then
									!print*,'cross-beta',i,j
                     							hplist(chnnum(i),chnnum(j)) = hplist(chnnum(i),chnnum(j))+1 
               							endif
                    					endif
                 				endif
              				endif

				enddo

			enddo
  
			ii = 1
			do i = 1,(nc+nc2)
				jj = 1
				ii = 1
				cluster(i,jj) = i
				layer(i,jj) = i
				do j = 1,(nc+nc2)
					if ((hblist(i,j).ge.3).or.(hplist(i,j).ge.3)) then
						jj = jj + 1
						cluster(i,jj) = j
					endif
					if (hblist(i,j).ge.3) then	
						ii = ii + 1
						layer(i,ii) = j
					endif
				enddo
			enddo
			count = 1
			do i = 1,(nc+nc2)
				do j = 1, size(pack(cluster(i,:), cluster(i,:)/= 0))
					if (link(cluster(i,j)) .eq. 0) then
						link(cluster(i,j)) = count
					else
						connect = link(cluster(i,j))
						do ii = 1, (nc+nc2)
							if (link(ii) .eq. connect) then
								link(ii) = count
							endif
						enddo
					endif
				enddo
			count = count + 1				
			enddo
			cluster = 0
			do i = 1, (nc+nc2)
				cluster(link(i),i) = i
			enddo
			linklayer = 0
			count = 1
			do i = 1,(nc+nc2)
				do j = 1, size(pack(layer(i,:), layer(i,:)/= 0))
					if (linklayer(layer(i,j)) .eq. 0) then
						linklayer(layer(i,j)) = count
					else
						connect = linklayer(layer(i,j))
						do ii = 1, (nc+nc2)
							if (linklayer(ii) .eq. connect) then
								linklayer(ii) = count
							endif
						enddo
					endif
				enddo
			count = count + 1				
			enddo
			
			layer = 0
			do i = 1, (nc+nc2)
				layer(linklayer(i),i) = i
			enddo

			do i = 1, (nc+nc2)
				clustersize = size(pack(cluster(i,:), cluster(i,:)/= 0))
				layersize = size(pack(layer(i,:), layer(i,:)/= 0))
				if (clustersize .eq. 1) then
					coil = coil + 1
				elseif ((clustersize .gt. 2).and.(layersize.lt.6)) then
					oligomer = oligomer + 1
					!write(8367,*) 'oligomer ', oligomer
					!write(8367,*) pack(cluster(i,:), cluster(i,:)/= 0)
				elseif ((clustersize .gt. 2).and.(layersize.ge.6)) then
					fibril = fibril + 1
					!write(8367,*) 'fibril ',fibril
					!write(8367,*) 'cluster'
					!write(8367,*) pack(cluster(i,:), cluster(i,:)/= 0)
					!write(8367,*) 'layer'
					!write(8367,*) pack(layer(i,:), layer(i,:)/= 0)
				endif	
			enddo
   	
end subroutine