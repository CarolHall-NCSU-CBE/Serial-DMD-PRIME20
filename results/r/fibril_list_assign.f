!find a list of peptides that are within the fibril
       subroutine fibril_list_assign()

#include "def.h"

        use global	
        implicit none

        integer i,j,ii,jj,a,aa,aaa,b,bb,bbb,c,nc
     	real*8 diff,corediscr,welldiscr,b_new,tij,depth
        real*8 dx,dy,dz,dvx,dvy,dvz,distance
		logical judge1,judge2
		real*8 dx1,dy1,dz1,dist1,dx2,dy2,dz2,dist2,angle
         
         nc = nop1/numbeads1+nop2/numbeads2    !Yiming
	     target_fibril = 0
         target_pep_num_in_fibril = 0	
         distance = 0.0		 
	
			  	do a=1,nc
            final_num_pep_in_fibril(a) = 0	
			final_num_sheet_in_fibril(a) = 0
            final_num_pep_in_sheet(a) = 0
				move(a) = 0
				move2(a)= 0
				sheet_identity(a) = 0
				pep_num(a) = 0
				pep_count(a) = 0
				pep_nbr_num(a) = 0
				num_sheet(a) = 0
				fibril_identity(a) = 0
				sheet_count(a) = 0
				sheet_nbr_num(a) = 0
				  do b=1,10
				pep_nbr(a,b) = 0
				sheet_nbr(a,b) = 0
				  enddo
			    do b=1,nc
			   temp_pep_in_sheet(a,b) = 0
			   temp_sheet_in_fibril(a,b) = 0
			   final_pep_in_fibril(a,b) = 0
			   final_pep_in_sheet(a,b) = 0
			   final_sheet_in_fibril(a,b) = 0
		        hb_contact(a,b) = 0
				hp_contact(a,b) = 0
				dimer(a,b)  = 0
				hp_dimer(a,b) = 0 	
				  contact_count(a,b) = 0				  
			   enddo 
			   enddo

		       do a =1,nc-1			   
		       do b = a+1,nc			   
                   do aa = (a-1)*numbeads1+chnln1+2,(a-1)*numbeads1+3*chnln1-1
				   do bb = (b-1)*numbeads1+chnln1+2,(b-1)*numbeads1+3*chnln1-1
	   	           if ((bptnr(aa).eq.bb).and.(chnnum(aa).ne.chnnum(bb))) then
			         hb_contact(a,b) = hb_contact(a,b) + 1	
			        !hb_contact(chnnum(bptnr(ii)),i) = hb_contact(chnnum(bptnr(ii)),i) + 1	
			       endif
		          enddo
                  enddo
					 do aa = (a-1)*numbeads1+3*chnln1+1,(a-1)*numbeads1+4*chnln1	 
					 do bb = (b-1)*numbeads1+3*chnln1+1,(b-1)*numbeads1+4*chnln1	
	         !dvx=sv(4,aa)-sv(4,bb)
	         !dvy=sv(5,aa)-sv(5,bb)	
	         !dvz=sv(6,aa)-sv(6,bb)
		    dx = rx(aa)-rx(bb)!+dvx*tfalse
		    dx = dx - dnint(dx)
            dy = ry(aa)-ry(bb)!+dvy*tfalse
		    dy = dy - dnint(dy)
            dz = rz(aa)-rz(bb)!+dvz*tfalse
	        dz = dz - dnint(dz)
            distance = sqrt(dx**2+dy**2+dz**2)
                if (distance*boxl.le.wel(identity(aa),identity(bb))) then 
			hp_contact(a,b) = hp_contact(a,b) + 1	
			hp_contact(b,a) = hp_contact(b,a) + 1	
			    endif				  				  
					 enddo
					 enddo				 
                if ((hp_contact(a,b).ge.chnln1/2).and.(hb_contact(a,b).eq.0)) then
			  hp_dimer(a,b) = 1
			  hp_dimer(b,a) = 1
			    endif
         enddo 
          enddo	
		  
		    do a =1,nc-1
		    do b = a+1,nc
			if (hb_contact(a,b).ge.chnln1/2+1) then
			dimer(a,b) = 1
			else if ((hb_contact(a,b).ge.1).and.previous_pep_in_fibril(a).and.previous_pep_in_fibril(b)) then
              dimer(a,b) = 1
            ! write(6,*) a,b			  
               endif			
			     enddo	
			   enddo   		

        !identify neighbour peptides 1 or 2 of each peptide
            do a=1,nc-1
               do b=a+1,nc
			     if (dimer(a,b).eq.1) then
			        pep_nbr_num(a) = pep_nbr_num(a) + 1
			        pep_nbr_num(b) = pep_nbr_num(b) + 1
					pep_nbr(a,pep_nbr_num(a)) = b
					pep_nbr(b,pep_nbr_num(b)) = a
			    endif
		       enddo
		    enddo
        !assign peptides to a sheet
			    sheet_num = 1		
		     do a=1,nc
		     if ((pep_nbr_num(a).ge.1).and.(sheet_identity(a).eq.0)) then	 ! doesn't matter 1 or 2
			    call sheet_assign(a)
			   do b=1,pep_num(sheet_num)
		      final_pep_in_sheet(sheet_num,b) = temp_pep_in_sheet(sheet_num,b)
                enddo	
           sheet_num = sheet_num + 1	! one sheet assignment finished          		
			   endif   
		     enddo
	 
                 sheet_num = sheet_num -1

!!!!!!! eliminate dimer and trimer from sheet !!!!!!!!!!!!!!!!!
				 
				 
              do a =1,nc
			  !final_sheet_identity(a) = sheet_identity(a) 
       if (sheet_identity(a).ne.0) final_num_pep_in_sheet(sheet_identity(a)) = final_num_pep_in_sheet(sheet_identity(a))+1					
		       enddo				 
				 
				  ! assigning sheets to fibril
           if(sheet_num.eq.1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            fibril_num = 1		  
			   fibril_identity(1) = fibril_num
			   num_sheet(fibril_num) = num_sheet(fibril_num) + 1    
			   final_sheet_in_fibril(fibril_num,num_sheet(fibril_num)) = 1 						
			 final_num_pep_in_fibril(fibril_identity(1)) = final_num_pep_in_fibril(fibril_identity(1))+pep_num(1)

				   do b = 1,nc
				    if (sheet_identity(b).ne.0) then 
					  if ((sheet_identity(b).ne.0).and.(pep_nbr_num(b).ge.1)) then					    
					    move(1) = move(1) + 1
				        final_pep_in_fibril(1,move(1)) = b
				      endif 	   
					 endif
                  enddo	
		
!!!!!!!!!  sheet_num = 1  !!!!!!!!!!!!!!		   
           elseif (sheet_num.gt.1) then 
				  do a = 1,sheet_num-1
				  do b = a+1,sheet_num			  
				     do aa = 1,pep_num(a)
				     do bb = 1,pep_num(b)
					 aaa = final_pep_in_sheet(a,aa)
					 bbb = final_pep_in_sheet(b,bb)
				  if (hp_dimer(aaa,bbb).eq.1) contact_count(a,b) = contact_count(a,b) + 1
					 enddo
					 enddo
				  
					 if (contact_count(a,b).ge.2)  then
                  sheet_nbr_num(a) = sheet_nbr_num(a) + 1
				  sheet_nbr_num(b) = sheet_nbr_num(b) + 1				  
				  sheet_nbr(a,sheet_nbr_num(a)) = b
				  sheet_nbr(b,sheet_nbr_num(b)) = a	
				  !write(6,*) i,j,contact_count(i,j),sheet_nbr_num(i),sheet_nbr_num(j)
				  endif
				  
					 enddo
                   enddo			   
				  
				   fibril_num = 1
				  
		     do a=1,sheet_num
		     if ((sheet_nbr_num(a).ne.0).and.(fibril_identity(a).eq.0)) then	 
			    call fibril_assign(a)
		         
                 fibril_num = fibril_num + 1		   
			   endif	
		     enddo
			 
				fibril_num = fibril_num -1
				
              do a =1,sheet_num 
			 final_num_pep_in_fibril(fibril_identity(a)) = final_num_pep_in_fibril(fibril_identity(a))+pep_num(a)		 
			   enddo		   
				   
		         a = 0
				   do b = 1,nc
				   a = fibril_identity(sheet_identity(b))
					  if ((sheet_identity(b).ne.0).and.(pep_nbr_num(b).ge.1)) then	 ! peptide at fibril end have HP contact with monomers
					    move2(a) = move2(a) + 1
				        final_pep_in_fibril(a,move2(a)) = b
				      endif 	   
                  enddo
	  
        endif  ! judge if sheet_num = 0,1,>2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do ii=1,fibril_num
        if (target_pep_num_in_fibril.le.final_num_pep_in_fibril(ii)) then 
		target_pep_num_in_fibril = final_num_pep_in_fibril(ii)
		target_fibril = ii
        endif
        enddo
		
		 do a=1,nc
			previous_pep_in_fibril(a) = .false.
		 enddo		
		
		if(sheet_num.eq.1) then
       do a=1,nc
       !    fibril_list(a) = final_pep_in_fibril(target_fibril,a)
        if (sheet_identity(a).eq.1)   previous_pep_in_fibril(a) = .true.
       enddo		
		elseif (sheet_num.gt.1) then
       do a=1,nc
       !    fibril_list(a) = final_pep_in_fibril(target_fibril,a)
	   if (target_fibril.eq.fibril_identity(sheet_identity(a))) previous_pep_in_fibril(a) = .true.
       enddo 
	   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       write(6,*) coll,pep_num(1),pep_num(2)
!       write(6,'(8i4)') fibril_num,sheet_num,num_sheet(1),final_num_pep_in_fibril(1),num_sheet(2),final_num_pep_in_fibril(2),num_sheet(3),final_num_pep_in_fibril(3)
	      
       end
	   
       recursive subroutine sheet_assign(x)
#include "def.h"	   
         use global
	      integer a,b,x,k,ii
 
			sheet_identity(x) = sheet_num
			pep_count(x) = pep_count(x) + 1

			   if (pep_count(x).eq.1)  then
			   pep_num(sheet_num) = pep_num(sheet_num) + 1
			temp_pep_in_sheet(sheet_num,pep_num(sheet_num)) = x
                   endif	
			  			   
			    if(pep_nbr_num(x).eq.1) then				
				  a = pep_nbr(x,1)

					 if(sheet_identity(a).eq.0) then
			sheet_identity(a) = sheet_num
			pep_count(a) = pep_count(a) + 1
			   if (pep_count(a).eq.1) then
			   pep_num(sheet_num) = pep_num(sheet_num) + 1
			temp_pep_in_sheet(sheet_num,pep_num(sheet_num)) = a 
              endif
			         endif

					  if(pep_nbr_num(a).gt.1) then
			!		     b = pep_nbr(a,2)
            !            if (sheet_identity(b).eq.0) call sheet_assign(b)
				     do ii=2,pep_nbr_num(a)
					  b = pep_nbr(a,ii)
			           if (sheet_identity(b).eq.0) call sheet_assign(b)
                      enddo				
					 endif
				  
	      elseif (pep_nbr_num(x).gt.1) then
				     do ii=1,pep_nbr_num(x)
					  a = pep_nbr(x,ii)
			           if (sheet_identity(a).eq.0) call sheet_assign(a)
                      enddo			   
                 endif			
         			
 		   return   

	     end
		 
		    recursive  subroutine fibril_assign(x) 
#include "def.h"	
           use global
	      integer a,b,x,k,ii

			fibril_identity(x) = fibril_num
			sheet_count(x) = sheet_count(x) + 1

			   if (sheet_count(x).eq.1) then
			   num_sheet(fibril_num) = num_sheet(fibril_num) + 1       ! be aware num_sheet() and sheet_num
			   temp_sheet_in_fibril(fibril_num,num_sheet(fibril_num)) = x 
              endif

			    if(sheet_nbr_num(x).eq.1) then				
				    a = sheet_nbr(x,1)

					 if(fibril_identity(a).eq.0) then
					      fibril_identity(a) = fibril_num
					      sheet_count(a) = sheet_count(a) + 1
					    if (sheet_count(a).eq.1) then
					     num_sheet(fibril_num) = num_sheet(fibril_num) + 1
					     temp_sheet_in_fibril(fibril_num,num_sheet(fibril_num)) = a 
				        endif
			         endif	  				 
                
				     if(sheet_nbr_num(a).gt.1) then			 
				!	      b = sheet_nbr(a,2)
                !       if (fibril_identity(b).eq.0) call fibril_assign(b)	
				     do ii=2,sheet_nbr_num(a)
					  b = sheet_nbr(a,ii)
			           if (fibril_identity(b).eq.0) call fibril_assign(b)
                      enddo					
					   endif
				  
				  elseif (sheet_nbr_num(x).gt.1) then
				     do ii=1,sheet_nbr_num(x)
					  a = sheet_nbr(x,ii)
			           if (fibril_identity(a).eq.0) call fibril_assign(a)
                      enddo		
				   		   
                 endif						  
			
          return
               
		     end

	

