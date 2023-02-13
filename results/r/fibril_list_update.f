!find a list of peptides that are within the fibril
       subroutine fibril_list_update()

#include "def.h"

        use global	
        implicit none

        integer ii,jj,a,aa,aaa,b,bb,bbb,c,nc,count,count2(nop1/numbeads1+nop2/numbeads2)
         logical judge
		 
         nc = nop1/numbeads1+nop2/numbeads2    !Yiming
			   
!!!!!!!!!!! global reassign dimer information
       if(target_pep_num_in_fibril.gt.0) then
               
			   do a=1,nc	
                count2(a) = 0
				move(a) = 0
				temp_fibril_list(a) = 0
				pep_nbr_num(a) = 0
				  do b=1,3
				pep_nbr(a,b) = 0
				  enddo
			    do b=1,nc
		        hb_contact(a,b) = 0
				dimer(a,b)  = 0				  
			    enddo 
			   enddo	   
			   
		       do a =1,nc-1			   
		       do b = a+1,nc			
                   do aa = (a-1)*numbeads1+chnln1+2,(a-1)*numbeads1+3*chnln1-1
				   do bb = (b-1)*numbeads1+chnln1+2,(b-1)*numbeads1+3*chnln1-1
	     if ((bptnr(aa).eq.bb).and.(chnnum(aa).ne.chnnum(bb)))  hb_contact(a,b) = hb_contact(a,b) + 1	
		          enddo
                  enddo
            enddo 
            enddo	
		  
		       do a =1,nc-1
		       do b = a+1,nc
              if (hb_contact(a,b).ge.chnln1/2+1)  dimer(a,b) = 1		  
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
		 
!!!!!!!!!!! 1. search within fibril list and delete
		 count = 0
         do a = 1,target_pep_num_in_fibril
		 b= fibril_list(a)
         if(pep_nbr_num(b).ge.1) then
          count = count + 1
		 temp_fibril_list(count) = b 
        endif
         enddo
		 target_pep_num_in_fibril = count
!!!!!!!!!!! 2. search monomer peptide list and add
         do a=1,nc
		 judge=.false.
           do b=1,count
		    c=temp_fibril_list(b)
         if(a.eq.c) judge=.true.
           enddo
		   
         if(judge) then
           else		
          do b=1,count
		    c=temp_fibril_list(b)		   
         if((pep_nbr(a,1).eq.c).or.(pep_nbr(a,2).eq.c).or.(pep_nbr(a,3).eq.c)) then
		   count2(a) = count2(a) + 1
         if (count2(a).eq.1) then
          target_pep_num_in_fibril  = count + 1
		   temp_fibril_list(target_pep_num_in_fibril) = a
	     endif
         endif
       enddo
	   
         endif
       enddo
		 
          do a = 1,nc
		 fibril_list(a)= 0
         enddo	      
	   
         do a = 1,target_pep_num_in_fibril
		 fibril_list(a)= temp_fibril_list(a)
         enddo		

	
       endif   	   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
!       write(6,*) coll,pep_num(1),pep_num(2)
!       write(6,'(8i4)') fibril_num,sheet_num,num_sheet(1),final_num_pep_in_fibril(1),num_sheet(2),final_num_pep_in_fibril(2),num_sheet(3),final_num_pep_in_fibril(3)
	      
       end