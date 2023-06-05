!	alex marchut
!	nc state university
!	feb. 7, 2003
!	this subroutine will check on the hydrogen bonding potential
!	if the particles are bound and the angle is bad then increment "boundbad"
!	if the particles are not bound and they should be- i.e. they are in the 
!	square well and the angle is good- then increment "unboundbad"

	subroutine check_nc_int(boundbad, unboundbad)

#include "def.h"

	use global
	use inputreadin

	implicit none

	integer evcode,type,i,j,ii,jj,boundbad,unboundbad,k,kstart,kend,m_ss,n_ss
	real*8 vxij,vyij,vzij,rxij,ryij,rzij,bij,rijsq,vijsq,d1,d2,d3,d4
	real*8 diff,corediscr,fac_sigsq,fac_cored,welldiscr,tij,rating

	n_ss=0
	m_ss=0

	!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do i=1,(noptotal)
		kstart=(i-1)*maxnbs+1
		kend=kstart+na_npt(i)-1
		do k=kstart,kend
              	j=nb(k)
	      		evcode=ev_code(j,i)
			if (evcode .eq. 15) then
                 		if (i .le. nop1) then
                 			ii=i-((chnnum(i)-1)*numbeads1)
                 		else
                 			ii=i-nop1-((chnnum(i)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
                 		endif
                 		if (j .le. nop1) then
                 			jj=j-((chnnum(j)-1)*numbeads1)
                 		else
                 			jj=j-nop1-((chnnum(j)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
                 		endif
        			if ((i .le. nop1) .and. (j .le. nop1)) then     
					if (bptnr(i) == j) then
!                   				i and j are already counted as forming a hydrogen bond
!		    				so check angle
!		    				if angle is bad something is fudged up
						if ((ii.ne.chnln1+1).and.(jj.ne.chnln1+1).and.(ii.ne.3*chnln1).and.(jj.ne.3*chnln1)) then
							if (identity(i)  .lt. identity(j)) then
!                         					i is n, j is c
		               				call repuls_check(i,j,rating)
                					else
!                         					i is c, j is n
			           				call repuls_check(j,i,rating)
                					endif
            					else
                       				rating=10.d0
            					endif
            					if (rating .gt. 10.d0) then
!		       				things are baaaaad!
		       				boundbad=boundbad+1
            					endif
						if ((ii.ne.chnln1+1).and.(jj.ne.chnln1+1).and.(ii.ne.3*chnln1).and.(jj.ne.3*chnln1)) then
		        				if (extra_repuls(i,4) .eq. j) then
		                   				n_ss = n_ss + 1
		            				else   
		                   				print*, 'no ss for bond', i , j
	                				endif		
	            				endif		
		 			else
		    				vxij=sv(4,i)-sv(4,j)
                    				vyij=sv(5,i)-sv(5,j)
                    				vzij=sv(6,i)-sv(6,j)
                    				rxij=sv(1,i)-sv(1,j)+vxij*tfalse
                    				ryij=sv(2,i)-sv(2,j)+vyij*tfalse
                    				rzij=sv(3,i)-sv(3,j)+vzij*tfalse
                    				rxij=rxij-dnint(rxij)
                    				ryij=ryij-dnint(ryij)
                    				rzij=rzij-dnint(rzij)
                    				bij=rxij*vxij+ryij*vyij+rzij*vzij
                    				rijsq=rxij*rxij+ryij*ryij+rzij*rzij
                    				vijsq=vxij*vxij+vyij*vyij+vzij*vzij
                    				diff=rijsq-welldia_sq(identity(i),identity(j))
                    				if (diff.lt.0.d0) then
!                   					then wells are already overlapping
!                      				so check angle
!                      				if angle is good something is fudged up
                          				if ((ii.ne.chnln1+1).and.(jj.ne.chnln1+1).and.(ii.ne.3*chnln1).and.(jj.ne.3*chnln1)) then
                          					if (identity(i) .lt. identity(j)) then
!                            					i is n, j is c
                             					call repuls_check(i,j,rating)
                          					else
!                            					i is c, j is n
                             					call repuls_check(j,i,rating)
                          					endif
                          				else
!                         					an end bead is involved, can't calc rating
!                         					dum...we'll let this slide
                          					rating = 11.d0
                          				endif
                          				if (rating.le.10.d0) then
!                         					the angle is good!
!                         					something is fudged up!
                          					unboundbad=unboundbad+1
                          				endif
		                  			if ((ii.ne.chnln1+1).and.(jj.ne.chnln1+1).and.(ii.ne.3*chnln1).and.(jj.ne.3*chnln1)) then
                          					if (extra_repuls(i,4) .eq. j) then
                             					n_ss = n_ss + 1
                          					else
                             					print*, 'no ss for non-bond', i , j
                          					endif
                          				endif
                      			endif
					endif
					if (extra_repuls(i,4) .eq. j) then
                    				m_ss = m_ss + 1
!                   				print*, coll, i, j 
                 			endif
				elseif ((i .le. nop1) .and. (j .gt. nop1)) then     
					if (bptnr(i) == j) then
!                   				i and j are already counted as forming a hydrogen bond
!		    				so check angle
!		    				if angle is bad something is fudged up
						if ((ii.ne.chnln1+1).and.(jj.ne.numbeads1+chnln2+1).and.(ii.ne.3*chnln1).and.(jj.ne.numbeads1+3*chnln2)) then
							if (identity(i)  .lt. identity(j)) then
!                         					i is n, j is c
		               				call repuls_check(i,j,rating)
                					else
!                         					i is c, j is n
			           				call repuls_check(j,i,rating)
                					endif
            					else
                       				rating=10.d0
            					endif
            					if (rating .gt. 10.d0) then
!		       				things are baaaaad!
		       				boundbad=boundbad+1
            					endif
						if ((ii.ne.chnln1+1).and.(jj.ne.numbeads1+chnln2+1).and.(ii.ne.3*chnln1).and.(jj.ne.numbeads1+3*chnln2)) then
		        				if (extra_repuls(i,4) .eq. j) then
		                   				n_ss = n_ss + 1
		            				else   
		                   				print*, 'no ss for bond', i , j
	                				endif		
	            				endif		
		 			else
		    				vxij=sv(4,i)-sv(4,j)
                    				vyij=sv(5,i)-sv(5,j)
                    				vzij=sv(6,i)-sv(6,j)
                    				rxij=sv(1,i)-sv(1,j)+vxij*tfalse
                    				ryij=sv(2,i)-sv(2,j)+vyij*tfalse
                    				rzij=sv(3,i)-sv(3,j)+vzij*tfalse
                    				rxij=rxij-dnint(rxij)
                    				ryij=ryij-dnint(ryij)
                    				rzij=rzij-dnint(rzij)
                    				bij=rxij*vxij+ryij*vyij+rzij*vzij
                    				rijsq=rxij*rxij+ryij*ryij+rzij*rzij
                    				vijsq=vxij*vxij+vyij*vyij+vzij*vzij
                    				diff=rijsq-welldia_sq(identity(i),identity(j))
                    				if (diff.lt.0.d0) then
!                   					then wells are already overlapping
!                      				so check angle
!                      				if angle is good something is fudged up
						if ((ii.ne.chnln1+1).and.(jj.ne.numbeads1+chnln2+1).and.(ii.ne.3*chnln1).and.(jj.ne.numbeads1+3*chnln2)) then
                          					if (identity(i) .lt. identity(j)) then
!                            					i is n, j is c
                             					call repuls_check(i,j,rating)
                          					else
!                            					i is c, j is n
                             					call repuls_check(j,i,rating)
                          					endif
                          				else
!                         					an end bead is involved, can't calc rating
!                         					dum...we'll let this slide
                          					rating = 11.d0
                          				endif
                          				if (rating.le.10.d0) then
!                         					the angle is good!
!                         					something is fudged up!
                          					unboundbad=unboundbad+1
                          				endif
						if ((ii.ne.chnln1+1).and.(jj.ne.numbeads1+chnln2+1).and.(ii.ne.3*chnln1).and.(jj.ne.numbeads1+3*chnln2)) then
                          					if (extra_repuls(i,4) .eq. j) then
                             					n_ss = n_ss + 1
                          					else
                             					print*, 'no ss for non-bond', i , j
                          					endif
                          				endif
                      			endif
					endif
					if (extra_repuls(i,4) .eq. j) then
                    				m_ss = m_ss + 1
!                   				print*, coll, i, j 
                 			endif                 
				elseif ((i .gt. nop1) .and. (j .le. nop1)) then     
					if (bptnr(i) == j) then
!                   				i and j are already counted as forming a hydrogen bond
!		    				so check angle
!		    				if angle is bad something is fudged up
						if ((ii.ne.numbeads1+chnln2+1).and.(jj.ne.chnln1+1).and.(ii.ne.numbeads1+3*chnln2).and.(jj.ne.3*chnln1)) then
							if (identity(i)  .lt. identity(j)) then
!                         					i is n, j is c
		               				call repuls_check(i,j,rating)
                					else
!                         					i is c, j is n
			           				call repuls_check(j,i,rating)
                					endif
            					else
                       				rating=10.d0
            					endif
            					if (rating .gt. 10.d0) then
!		       				things are baaaaad!
		       				boundbad=boundbad+1
            					endif
						if ((ii.ne.numbeads1+chnln2+1).and.(jj.ne.chnln1+1).and.(ii.ne.numbeads1+3*chnln2).and.(jj.ne.3*chnln1)) then
		        				if (extra_repuls(i,4) .eq. j) then
		                   				n_ss = n_ss + 1
		            				else   
		                   				print*, 'no ss for bond', i , j
	                				endif		
	            				endif		
		 			else
		    				vxij=sv(4,i)-sv(4,j)
                    				vyij=sv(5,i)-sv(5,j)
                    				vzij=sv(6,i)-sv(6,j)
                    				rxij=sv(1,i)-sv(1,j)+vxij*tfalse
                    				ryij=sv(2,i)-sv(2,j)+vyij*tfalse
                    				rzij=sv(3,i)-sv(3,j)+vzij*tfalse
                    				rxij=rxij-dnint(rxij)
                    				ryij=ryij-dnint(ryij)
                    				rzij=rzij-dnint(rzij)
                    				bij=rxij*vxij+ryij*vyij+rzij*vzij
                    				rijsq=rxij*rxij+ryij*ryij+rzij*rzij
                    				vijsq=vxij*vxij+vyij*vyij+vzij*vzij
                    				diff=rijsq-welldia_sq(identity(i),identity(j))
                    				if (diff.lt.0.d0) then
!                   					then wells are already overlapping
!                      				so check angle
!                      				if angle is good something is fudged up
							if ((ii.ne.numbeads1+chnln2+1).and.(jj.ne.chnln1+1).and.(ii.ne.numbeads1+3*chnln2).and.(jj.ne.3*chnln1)) then
                          					if (identity(i) .lt. identity(j)) then
!                            					i is n, j is c
                             					call repuls_check(i,j,rating)
                          					else
!                            					i is c, j is n
                             					call repuls_check(j,i,rating)
                          					endif
                          				else
!                         					an end bead is involved, can't calc rating
!                         					dum...we'll let this slide
                          					rating = 11.d0
                          				endif
                          				if (rating.le.10.d0) then
!                         					the angle is good!
!                         					something is fudged up!
                          					unboundbad=unboundbad+1
                          				endif
							if ((ii.ne.numbeads1+chnln2+1).and.(jj.ne.chnln1+1).and.(ii.ne.numbeads1+3*chnln2).and.(jj.ne.3*chnln1)) then
                          					if (extra_repuls(i,4) .eq. j) then
                             					n_ss = n_ss + 1
                          					else
                             					print*, 'no ss for non-bond', i , j
                          					endif
                          				endif
                      			endif
					endif
					if (extra_repuls(i,4) .eq. j) then
                    				m_ss = m_ss + 1
!                   				print*, coll, i, j 
                 			endif                 
				elseif ((i .gt. nop1) .and. (j .gt. nop1)) then     
					if (bptnr(i) == j) then
!                   				i and j are already counted as forming a hydrogen bond
!		    				so check angle
!		    				if angle is bad something is fudged up
						if ((ii.ne.numbeads1+chnln2+1).and.(jj.ne.numbeads1+chnln2+1).and.(ii.ne.numbeads1+3*chnln2).and.(jj.ne.numbeads1+3*chnln2)) then
							if (identity(i)  .lt. identity(j)) then
!                         					i is n, j is c
		               				call repuls_check(i,j,rating)
                					else
!                         					i is c, j is n
			           				call repuls_check(j,i,rating)
                					endif
            					else
                       				rating=10.d0
            					endif
            					if (rating .gt. 10.d0) then
!		       				things are baaaaad!
		       				boundbad=boundbad+1
            					endif
						if ((ii.ne.numbeads1+chnln2+1).and.(jj.ne.numbeads1+chnln2+1).and.(ii.ne.numbeads1+3*chnln2).and.(jj.ne.numbeads1+3*chnln2)) then
		        				if (extra_repuls(i,4) .eq. j) then
		                   				n_ss = n_ss + 1
		            				else   
		                   				print*, 'no ss for bond', i , j
	                				endif		
	            				endif		
		 			else
		    				vxij=sv(4,i)-sv(4,j)
                    				vyij=sv(5,i)-sv(5,j)
                    				vzij=sv(6,i)-sv(6,j)
                    				rxij=sv(1,i)-sv(1,j)+vxij*tfalse
                    				ryij=sv(2,i)-sv(2,j)+vyij*tfalse
                    				rzij=sv(3,i)-sv(3,j)+vzij*tfalse
                    				rxij=rxij-dnint(rxij)
                    				ryij=ryij-dnint(ryij)
                    				rzij=rzij-dnint(rzij)
                    				bij=rxij*vxij+ryij*vyij+rzij*vzij
                    				rijsq=rxij*rxij+ryij*ryij+rzij*rzij
                    				vijsq=vxij*vxij+vyij*vyij+vzij*vzij
                    				diff=rijsq-welldia_sq(identity(i),identity(j))
                    				if (diff.lt.0.d0) then
!                   					then wells are already overlapping
!                      				so check angle
!                      				if angle is good something is fudged up
							if ((ii.ne.numbeads1+chnln2+1).and.(jj.ne.numbeads1+chnln2+1).and.(ii.ne.numbeads1+3*chnln2).and.(jj.ne.numbeads1+3*chnln2)) then
                          					if (identity(i) .lt. identity(j)) then
!                            					i is n, j is c
                             					call repuls_check(i,j,rating)
                          					else
!                            					i is c, j is n
                             					call repuls_check(j,i,rating)
                          					endif
                          				else
!                         					an end bead is involved, can't calc rating
!                         					dum...we'll let this slide
                          					rating = 11.d0
                          				endif
                          				if (rating.le.10.d0) then
!                         					the angle is good!
!                         					something is fudged up!
                          					unboundbad=unboundbad+1
                          				endif
							if ((ii.ne.numbeads1+chnln2+1).and.(jj.ne.numbeads1+chnln2+1).and.(ii.ne.numbeads1+3*chnln2).and.(jj.ne.numbeads1+3*chnln2)) then
                          					if (extra_repuls(i,4) .eq. j) then
                             					n_ss = n_ss + 1
                          					else
                             					print*, 'no ss for non-bond', i , j
                          					endif
                          				endif
                      			endif
					endif
					if (extra_repuls(i,4) .eq. j) then
                    				m_ss = m_ss + 1
!                   				print*, coll, i, j 
                 			endif
				endif
			endif                 
	   	enddo
	enddo

	if (m_ss .ne. n_ss) then
		print*, 'wrong ss at coll', coll, n_ss, m_ss
		call exit(-1)
	endif

	return

	end
	

