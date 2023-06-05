!	subroutine make_code generates the matrix ev_code(nop,nop) 
!	that is used to for looking up identities (row 1) and event types
!	based on i,j pair also generates ev_param file (row 1: "factor", 
!	row 2: "blmin", and row 3: "blmax")

	subroutine make_code()
 
#include "def.h"

	use global
	use inputreadin

	implicit none

	integer i,j,k,kk,l,ll,lll,m,mm,ncount

!	make parameter matrix

	do k=1,3
		do l=1,50
			ev_param(k,l)=1
		enddo
	enddo

!	ev_param(1,15):  so oh in bond don't pass over each other
!	2.24 comes from length of co bond + length of nh bond
!	1.05 is just to make sure core collision happens before o and h
!	get so close that they can cross paths
	ev_param(1,15)=1.05d0*((2.24d0/boxl_orig)/((sigma(1)+sigma(4))/2.d0))
	ev_param(1,17)=sqz1
	ev_param(1,18)=sqz2
	ev_param(1,19)=sqz3
	ev_param(1,20)=sqz4
	ev_param(1,21)=sqz5
!	ev_param(1,22)=sqz6
!	ev_param(1,23)=sqz7
!	ev_param(1,24)=sqz8
!	ev_param(1,25)=sqz9
!	ev_param(1,26)=sqz10
	ev_param(1,27)=sqz11
!	ev_code = 4:  dnc = c_alpha_i to n_i (covalent bond)
!	ev_code = 5:  dcc = c_alpha_i to c_i (covalent bond)
!	ev_code = 6:  dcn = c_i to n_i+1, peptide bond (covalent bond)
!   	ev_code = 7:  dtie = c_alpha_i to n_i+1 (pseudobond)
!   	ev_code = 8:  dtie2 = c_i to c_alpha_i+1 and c_i to n_i (pseudobond)
!   	ev_code = 9:  dcaca = c_alpha_i to c_alpha_i+1 (pseudobond)
!   	ev_code = 10: c_alpha_i to r_i (covalent bond)
!     	ev_code = 11: r_i to n_i (pseudobond)
!    	ev_code = 12: r_i to c_i (pseudobond) 
  	ev_param(2,4)=dnc*(1.d0-del)/boxl_orig	!ev_code = 4
     	ev_param(2,5)=dcc*(1.d0-del)/boxl_orig	!ev_code = 5
 	ev_param(2,6)=dcn*(1.d0-del)/boxl_orig	! ev_code = 6
  	ev_param(2,7)=dtie*(1.d0-del)/boxl_orig	! ev_code = 7
  	ev_param(2,8)=dtie2*(1.d0-del)/boxl_orig	! ev_code = 8
    	ev_param(2,9)=dcaca*(1.d0-del)/boxl_orig	! ev_code = 9
    	ev_param(2,10)=(1.d0-del)	
     	ev_param(2,11)=(1.d0-del)
    	ev_param(2,12)=(1.d0-del)
   	ev_param(3,4)=dnc*(1.d0+del)/boxl_orig
   	ev_param(3,5)=dcc*(1.d0+del)/boxl_orig
   	ev_param(3,6)=dcn*(1.d0+del)/boxl_orig
   	ev_param(3,7)=(1.d0+del)*dtie/boxl_orig
   	ev_param(3,8)=(1.d0+del)*dtie2/boxl_orig
    	ev_param(3,9)=(1.d0+del)*dcaca/boxl_orig
    	ev_param(3,10)=(1.d0+del)
   	ev_param(3,11)=(1.d0+del)
     	ev_param(3,12)=(1.d0+del)
!     	note:  for ev_code=10 (bond event between side chain and its
!   	c alpha, must later multiply these blmin and blmax values by the 
!    	length of that particular bond
!  	default event type is core collision

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do k=1,(noptotal)
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
		do l=1,(noptotal)
			ev_code(k,l)=1
		enddo
	enddo
        
!  	interactions between all r's - set to 16 if both hp()=1

!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do m=1,(noptotal)-1
		if (m .le. nop1) then
			mm=m-(chnnum(m)-1)*numbeads1
!LR: Changed an open else statement to a constrained else-if statement (m in species 2)
		elseif (m .le. nop1+nop2) then
          		mm=m-nop1-((chnnum(m)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
!LR: Created a constrained else-if statement (m in species 3) 

          	endif
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
           	do l=m+1,(noptotal)
	    		if (l .le. nop1) then
	    			ll=l-(chnnum(l)-1)*numbeads1
!LR: Changed an open else statement to a constrained else-if statement (m in species 2)
           		elseif (l .le. nop1+nop2) then
           			ll=l-nop1-((chnnum(l)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
!LR: Created a constrained else-if statement (m in species 3) 

           		endif
              	if ((identity(m).gt.8).and.(identity(l).gt.8)) then
                 		if ((hp(mm).eq.1).and.(hp(ll).eq.1)) then
                    			if ((chnnum(l).ne.chnnum(m))) then
                       			ev_code(m,l)=16
                    			endif
                 		endif
              	endif
           	enddo
	enddo

!    	interactions between all nc's - set to 15 

#ifndef no_hbs
	
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
	do m=1,(noptotal)-1
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
		do l=m+1,(noptotal)
              	if (((identity(m).eq.1).and.(identity(l).eq.4)).or. &
&                 	((identity(m).eq.4).and.(identity(l).eq.1))) then
                    		ev_code(m,l)=15
              	endif
           	enddo
	enddo

!   	reset same-chain nc's and hp's to core
	do k=1,numchains-(nop2/numbeads2)
		kk=(k-1)*numbeads1
           	do m=kk+1,kk+numbeads1-1
              	do l=m+1,kk+numbeads1
                 		ev_code(m,l)=1
              	enddo
           	enddo
	enddo
	do k=(nop1/numbeads1)+1,numchains
		kk= nop1+(k-(nop1/numbeads1)-1)*numbeads2
		do m=kk+1,kk+numbeads2-1
			do l=m+1,kk+numbeads2
				ev_code(m,l)=1
			enddo
		enddo
	enddo
	
! Make sure proline NH doesn't HB
        DO m=1,noptotal-1
           DO l=m+1,noptotal	
        		if ((m .le. nop1) .and. (l .le. nop1)) then
					IF ((identity(m).eq.17).and.(identity(l).eq.4)) THEN
						ev_code(m-chnln1*2,l)=1
					ELSE IF ((identity(m).eq.4).and.(identity(l).eq.17)) THEN
						ev_code(l-chnln1*2,m)=1
					ENDIF
      			elseif ((m .le. nop1) .and. (l .gt. nop1)) then
					IF ((identity(m).eq.17).and.(identity(l).eq.4)) THEN
						ev_code(m-chnln1*2,l)=1
					ELSE IF ((identity(m).eq.4).and.(identity(l).eq.17)) THEN
						ev_code(l-chnln2*2,m)=1
					ENDIF
                elseif ((m .gt. nop1) .and. (l .gt. nop1)) then
					IF ((identity(m).eq.17).and.(identity(l).eq.4)) THEN
						ev_code(m-chnln2*2,l)=1
					ELSE IF ((identity(m).eq.4).and.(identity(l).eq.17)) THEN
						ev_code(l-chnln2*2,m)=1
					ENDIF
                endif              
           ENDDO
        ENDDO	
	
#else
        write(fileout,*)'****   no hydrogen bonding in this run ****'
#endif

	do ll=1,(nop1/numbeads1)
		lll=(ll-1)*numbeads1
!          	initializing event types these are ints btwn all r pairs 
!          	except pairs involving first and/or last r
           	ncount=2+n_b_hydro
           	do k=1,chnln1
              	if (fside1(k) .ne. 0) then
                 		do l=ncount,chnln1
					if (fside1(l) .ne. 0) then
                       			if ((hp1(fside1(k)).eq.1).and.(hp1(fside1(l)).eq.1)) then
                          				ev_code(lll+fside1(k),lll+fside1(l))=16
                       			endif
                    			endif
                 		enddo
              	endif
			ncount=ncount+1
		enddo

#ifndef no_hbs
!   		set the hbonding pairs
           	do k=lll+chnln1+1,lll+2*chnln1
              	do l=lll+2*chnln1+1,lll+3*chnln1
                 		ev_code(k,l)=15
              	enddo
           	enddo

!     		undo sqwel for n,c, too close

		do k=lll+chnln1+1,lll+2*chnln1
              	do l=k+chnln1-n_b_hbond,k+chnln1+n_b_hbond
                 		ev_code(k,l)=1
              	enddo
           	enddo
#endif

!          	c_alpha_i to n_i

           	do k=lll+1,lll+chnln1
              	ev_code(k,chnln1+k)=4
           	enddo

!          	c_alpha_i to c_i

           	do k=lll+1,lll+chnln1
              	ev_code(k,2*chnln1+k)=5 
           	enddo

!         	c_i to n_i+1, peptide bond 
 
           	do k=lll+chnln1+2,lll+2*chnln1
             		ev_code(k,k+chnln1-1)=6
           	enddo

!          	c_alpha_i to n_i+1

           	do k=lll+1,lll+chnln1-1
              	ev_code(k,chnln1+k+1)=7
           	enddo
           
!          	c_i to c_alpha_i+1  
  
           	do k=lll+2,lll+chnln1
              	ev_code(k,2*chnln1+k-1)=8
           	enddo

!          	c_i to n_i

           	do k=lll+chnln1+1,lll+2*chnln1
              	ev_code(k,k+chnln1)=8
           	enddo

!          	c_alpha_i to c_alpha_i+1

           	do k=lll+1,lll+chnln1-1
              	ev_code(k,k+1)=9
           	enddo

!          	c_alpha_i to r_i

           	do k=lll+1,lll+chnln1
              	if (fside1(k-lll) .ne. 0) ev_code(k,lll+fside1(k-lll))=10
           	enddo

!          	connect r_i to n_i = bond drn

           	do k=lll+chnln1+1,lll+2*chnln1
              	if (fside1(k-lll-chnln1) .ne. 0) ev_code(k,lll+fside1(k-lll-chnln1))=11
           	enddo

!          	connect r_i to c_i = bond drc

           	do k=lll+2*chnln1+1,lll+3*chnln1
              	if (fside1(k-lll-2*chnln1) .ne. 0) ev_code(k,lll+fside1(k-lll-2*chnln1))=12
           	enddo

!          	c_alpha_i to c_i+1,
 
           	do k=lll+1,lll+chnln1-1
              	ev_code(k,k+2*chnln1+1)=17
           	enddo

!          	c_alpha_i to n_i-1,
 
           	do k=lll+2,lll+chnln1
              	ev_code(k,k+chnln1-1)=18
           	enddo

!          	c_i to n_i+2

           	do k=lll+chnln1+3,lll+2*chnln1
              	ev_code(k,k+chnln1-2)=19
           	enddo

!          	n_i to n_i+1

           	do k=lll+chnln1+1,lll+2*chnln1-1
              	ev_code(k,k+1)=20
           	enddo

!          	c_i to c_i+1

           	do k=lll+2*chnln1+1,lll+3*chnln1-1
              	ev_code(k,k+1)=21
           	enddo
        
!          	r_i to c_i-1

           	do k=lll+2*chnln1+1,lll+3*chnln1-1
              	if (fside1(k-lll-2*chnln1+1) .ne. 0) ev_code(k,lll+fside1(k-lll-2*chnln1+1))=22
           	enddo

!          	r_i to n_i+1

           	do k=lll+chnln1+2,lll+2*chnln1
              	if (fside1(k-lll-chnln1-1) .ne. 0) ev_code(k,lll+fside1(k-lll-chnln1-1))=23
        	enddo

!          	r_i to c_alpha_i-1

           	do k=lll+1,lll+chnln1-1
              	if (fside1(k-lll+1) .ne. 0) ev_code(k,lll+fside1(k-lll+1))=24
           	enddo

!          	r_i to c_alpha_i+1

           	do k=lll+2,lll+chnln1
              	if (fside1(k-lll-1) .ne. 0) ev_code(k,lll+fside1(k-lll-1))=25
           	enddo

!          	c_i-1 to r_i+1

           	do k=lll+2*chnln1+1,lll+3*chnln1-2
              	if (fside1(k-lll-2*chnln1+2) .ne. 0) ev_code(k,lll+fside1(k-lll-2*chnln1+2))=26
           	enddo
	enddo

	do ll=1,(nop2/numbeads2)
		lll=nop1+(ll-1)*numbeads2
!          	initializing event types these are ints btwn all r pairs 
!          	except pairs involving first and/or last r
           	ncount=2+n_b_hydro
           	do k=1,chnln2
              	if (fside2(k) .ne. 0) then
                 		do l=ncount,chnln2
                    			if (fside2(l) .ne. 0) then
                       			if ((hp2(fside2(k)-nop1).eq.1).and.(hp2(fside2(l)-nop1).eq.1)) then
                          				ev_code(lll+fside2(k)-nop1,lll+fside2(l)-nop1)=16
                       			endif
                    			endif
                 		enddo
              	endif
              	ncount=ncount+1
           	enddo
        
#ifndef no_hbs
	if(chaptype .eq. 1) then
!          	set the hbonding pairs
           	do k=lll+chnln2+1,lll+2*chnln2
              	do l=lll+2*chnln2+1,lll+3*chnln2
                 		ev_code(k,l)=15
              	enddo
           	enddo

!          	undo sqwel for n,c, too close

           	do k=lll+chnln2+1,lll+2*chnln2
              	do l=k+chnln2-n_b_hbond,k+chnln2+n_b_hbond
                 		ev_code(k,l)=1
              	enddo
           	enddo
	endif
#endif
	if(chaptype .eq. 1) then
!          	c_alpha_i to n_i

           	do k=lll+1,lll+chnln2
              	ev_code(k,chnln2+k)=4
           	enddo

!          	c_alpha_i to c_i

           	do k=lll+1,lll+chnln2
              	ev_code(k,2*chnln2+k)=5 
           	enddo

!          	c_i to n_i+1, peptide bond   

           	do k=lll+chnln2+2,lll+2*chnln2
            		ev_code(k,k+chnln2-1)=6
           	enddo

!          	c_alpha_i to n_i+1

           	do k=lll+1,lll+chnln2-1
              	ev_code(k,chnln2+k+1)=7
           	enddo
           
!          	c_i to c_alpha_i+1  
  
           	do k=lll+2,lll+chnln2
              	ev_code(k,2*chnln2+k-1)=8
           	enddo

!          	c_i to n_i

           	do k=lll+chnln2+1,lll+2*chnln2
              	ev_code(k,k+chnln2)=8
           	enddo

!          	c_alpha_i to c_alpha_i+1

           	do k=lll+1,lll+chnln2-1
              	ev_code(k,k+1)=9
           	enddo

!          	c_alpha_i to r_i

           	do k=lll+1,lll+chnln2
              	if (fside2(k-lll) .ne. 0) ev_code(k,lll+fside2(k-lll)-nop1)=10
           	enddo

!          	connect r_i to n_i = bond drn

           	do k=lll+chnln2+1,lll+2*chnln2
              	if (fside2(k-lll-chnln2) .ne. 0) ev_code(k,lll+fside2(k-lll-chnln2)-nop1)=11
           	enddo

!          	connect r_i to c_i = bond drc

           	do k=lll+2*chnln2+1,lll+3*chnln2
              	if (fside2(k-lll-2*chnln2) .ne. 0) ev_code(k,lll+fside2(k-lll-2*chnln2)-nop1)=12
           	enddo

!          	c_alpha_i to c_i+1, 

           	do k=lll+1,lll+chnln2-1
              	ev_code(k,k+2*chnln2+1)=17
           	enddo

!          	c_alpha_i to n_i-1, 

           	do k=lll+2,lll+chnln2
              	ev_code(k,k+chnln2-1)=18
           	enddo

!          	c_i to n_i+2

           	do k=lll+chnln2+3,lll+2*chnln2
              	ev_code(k,k+chnln2-2)=19
           	enddo

!          	n_i to n_i+1

           	do k=lll+chnln2+1,lll+2*chnln2-1
              	ev_code(k,k+1)=20
           	enddo

!          	c_i to c_i+1

           	do k=lll+2*chnln2+1,lll+3*chnln2-1
              	ev_code(k,k+1)=21
           	enddo
        
!          	r_i to c_i-1
	
           	do k=lll+2*chnln2+1,lll+3*chnln2-1
              	if (fside2(k-lll-2*chnln2+1) .ne. 0) ev_code(k,lll+fside2(k-lll-2*chnln2+1)-nop1)=22
           	enddo

!          	r_i to n_i+1

           	do k=lll+chnln2+2,lll+2*chnln2
              	if (fside2(k-lll-chnln2-1) .ne. 0) ev_code(k,lll+fside2(k-lll-chnln2-1)-nop1)=23
           	enddo

!          	r_i to c_alpha_i-1

           	do k=lll+1,lll+chnln2-1
              	if (fside2(k-lll+1) .ne. 0) ev_code(k,lll+fside2(k-lll+1)-nop1)=24
           	enddo

!          	r_i to c_alpha_i+1

           	do k=lll+2,lll+chnln2
              	if (fside2(k-lll-1) .ne. 0) ev_code(k,lll+fside2(k-lll-1)-nop1)=25
           	enddo

!          	c_i-1 to r_i+1

           	do k=lll+2*chnln2+1,lll+3*chnln2-2
              	if (fside2(k-lll-2*chnln2+2) .ne. 0) ev_code(k,lll+fside2(k-lll-2*chnln2+2)-nop1)=26
           	enddo
	endif
	enddo

	if (chaptype .eq. 2) then
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
		do m=1,(noptotal)-1
			if (m .le. nop1) then
				mm=m-(chnnum(m)-1)*numbeads1
!LR: Changed an open else statement to a constrained else-if statement (m in species 2)
			elseif (m .le. nop1+nop2) then
          			mm=m-nop1-((chnnum(m)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
!LR: Created a constrained else-if statement (m in species 3) 
          		endif
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
           		do l=m+1,(noptotal)
	    			if (l .le. nop1) then
	    				ll=l-(chnnum(l)-1)*numbeads1
!LR: Changed an open else statement to a constrained else-if statement (m in species 2)
           			elseif (l .le. nop1+nop2) then
           				ll=l-nop1-((chnnum(l)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
!LR: Created a constrained else-if statement (m in species 3) 
           			endif
              		if ((identity(m).gt.8).and.(identity(l).gt.8)) then
                 			if ((hp(mm).eq.1).and.(hp(ll).eq.1)) then
                    				if ((chnnum(l).ne.chnnum(m))) then
                       				ev_code(m,l)=16
                    				endif
                 			endif
              		endif
           		enddo
		enddo	
	endif	

	if (chaptype .eq. 2) then
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
		do m=nop1+1,(noptotal)-1
			if (m .le. nop1) then
				mm=m-(chnnum(m)-1)*numbeads1
!LR: Changed an open else statement to a constrained else-if statement (m in species 2)
			elseif (m .le. nop1+nop2) then
          			mm=m-nop1-((chnnum(m)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
!LR: Created a constrained else-if statement (m in species 3) 
          		endif
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
           		do l=m+1,(noptotal)
	    			if (l .le. nop1) then
	    				ll=l-(chnnum(l)-1)*numbeads1
!LR: Changed an open else statement to a constrained else-if statement (m in species 2)
           			elseif (l .le. nop1+nop2) then
           				ll=l-nop1-((chnnum(l)-(nop1/numbeads1)-1)*numbeads2)+numbeads1
!LR: Created a constrained else-if statement (m in species 3)            			
           			endif
              		if ((identity(m).gt.8).and.(identity(l).gt.8)) then
                 			if ((hp(mm).eq.1).and.(hp(ll).eq.1)) then
                    				if ((chnnum(l).ne.chnnum(m))) then
                       				ev_code(m,l)=1
                    				endif
                 			endif
              		endif
           		enddo
		enddo	
	endif	


!LR: This comment applies to all code between LRStart and LREnd - 
!LR: This block is an exact duplicate of the 1st and 2nd species make code
! block. Again, I exactly copied the  2nd species block (which was added by Dave Latshaw), 
! I just changed the variable bounds where necessary.
!LRStart - start of 3rd species block
	
!LREnd - end of 3rd species block
	
	
!  	put ev_code_orig in unused triangular section of ev_code
!LR: Changed a hardcoded 2-species variable reference to a noptotal variable
   	do i=1,(noptotal)-1
   !LR: Changed a hardcoded 2-species variable reference to a noptotal variable
    		do j=i+1,(noptotal)
              	ev_code(j,i)=ev_code(i,j)
           	enddo
   	enddo

   	return

   	end
