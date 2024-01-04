	subroutine nbor_setup()

#include "def.h"

	use global
	use inputreadin

	implicit none

	real*8 sig_ij,sig_max(50),sig,sigij,sigi,sigj,welli,wellj,rl
	integer evcode,i,j

!	set up nbor list as function of bead size
	sig_max_all = 0.d0

      	do i = 1, 50
		sig_max(i) = 0.d0
      	enddo

		!LR: species 1 vs species 1
      	do i = 1, numbeads1 - 1
		do j = i + 1, numbeads1
            		welli=welldia(identity(i))
            		wellj=welldia(identity(j))
	    		evcode=ev_code(j,i)
	            	if (evcode .le. 26) then 
               		if (evcode.le.3) then
	          			sig_ij=sigma_2b(identity(i),identity(j))
               		elseif (evcode.eq.15) then
	          			sig_ij=0.5d0*(welli+wellj)
	       		elseif (evcode.eq.16) then
                  			sig_ij=wel(identity(i),identity(j))/boxl_orig
	       		else
	          			sig_ij=0.d0
               		endif
	       		if (sig_ij .gt. sig_max(evcode)) then
                  			sig_max(evcode) = sig_ij
	          			if (sig_max(evcode) .gt. sig_max_all) sig_max_all = sig_max(evcode)
               		endif
	    		end if
         	end do
	end do

	!LR: species 2 vs. species 2
	do i = nop1+1, nop1+numbeads2 - 1
		do j = i + 1, nop1+numbeads2
            		welli=welldia(identity(i))
            		wellj=welldia(identity(j))
	    		evcode=ev_code(j,i)
            		if (evcode .le. 26) then 
               		if (evcode.le.3) then
	          			sig_ij=sigma_2b(identity(i),identity(j))
               		elseif (evcode.eq.15) then
	          			sig_ij=0.5d0*(welli+wellj)
	       		elseif (evcode.eq.16) then
                  			sig_ij=wel(identity(i),identity(j))/boxl_orig 
	       		else
	          			sig_ij=0.d0
               		endif
	       		if (sig_ij .gt. sig_max(evcode)) then
                  			sig_max(evcode) = sig_ij
	          			if (sig_max(evcode) .gt. sig_max_all) sig_max_all = sig_max(evcode)
               		endif
	    		end if
         	end do
	end do

	!LR: species 1 vs. species 2
	do i = 1, numbeads1
		do j = nop1+1, nop1+numbeads2
            		welli=welldia(identity(i))
            		wellj=welldia(identity(j))
	   		evcode=ev_code(j,i)
            		if (evcode .le. 26) then 
               		if (evcode.le.3) then
	          			sig_ij=sigma_2b(identity(i),identity(j))
               			elseif (evcode.eq.15) then
	          			sig_ij=0.5d0*(welli+wellj)
	       		elseif (evcode.eq.16) then
                  			sig_ij=wel(identity(i),identity(j))/boxl_orig 
	       		else
	          			sig_ij=0.d0
               		endif
	       		if (sig_ij .gt. sig_max(evcode)) then
                  			sig_max(evcode) = sig_ij
	          			if (sig_max(evcode) .gt. sig_max_all) sig_max_all = sig_max(evcode)
               		endif
	    		end if
         	end do
	end do
	
	!LR: species 1 vs. species 3 - functionally identical to Dave Latshaw's species 1 vs. species 2 code.

	!species 2 vs. species 3 - functionally identical to Dave Latshaw's species 1 vs. species 2 code

	!LR: no species 3 vs. species 3 b/c there will only ever be one NP

	sigij=0.d0

      	do i = 1, 4
		do j = 1, 4
	    		if ((i .ne. 3) .or. (j .ne. 3)) sig=sigma(i)+sigma(j)
	    		if (sig .gt. sigij) sigij = sig
         	end do
	enddo

      	do i = 40, 50
	 	sig_ij=0.5d0*sigij*1.d0
         	sig_max(i) = sig_ij
	 	if (sig_max(i) .gt. sig_max_all) sig_max_all = sig_max(i)
      	enddo

      	do i = 1, 50
	 	rlsq(i) = 0.d0
	 	if (sig_max(i) .ne. 0.d0) rlsq(i) = ((rl_const-1)*sig_max_all+sig_max(i))**2
      	enddo

      	rl=rl_const*sig_max_all
      	hdelr=(0.4d0*(rl-sig_max_all))**2

#ifdef glycine 
      	do i = 1, nop
         	bdln_dummy(i)=bdln(i-numbeads*(chnnum(i)-1))
      	end do
#endif

      	write(6,*)' '
76    	format(' rl_const = ',f8.4,'  sig_max_all = ',f8.4)
      	write(6,76)rl_const,sig_max_all*boxl_orig

      	return

      	end
