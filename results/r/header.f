      module global
#include "def.h"	
       real n_phi_shl(2),n_psi_shl(2),pi
	integer*8 coll,ncoll	
        real*8 t,ttotal,t_old,t_temp
        integer ab(numbeads1+numbeads2),identity(nop1+nop2),fside1(chnln1),residue(nop1+nop2),chnnum(nop1+nop2),fside2(chnln2)
        real*8  eptemp,ep(28,28),ep2(28,28),bdtemp,wltemp,wel(28,28),wel2(28,28),bds(28,28)
        real*8 rx(nop1+nop2),ry(nop1+nop2),rz(nop1+nop2)
        integer noptotal,bptnr(nop1+nop2),count_file
		real*8 hb_total,std1,hb2_total,std2,hb_cross_total,std3
      	real*8 ered,tred,rg_avg,e2e_avg,ehh_ii,ehh_ij
      	integer hb_alpha,hb_ii,hb_ij		

! Yiming		
		integer contact_count(nop1/numbeads1+nop2/numbeads2,nop1/numbeads1+nop2/numbeads2),hb_contact(nop1/numbeads1+nop2/numbeads2,nop1/numbeads1+nop2/numbeads2),hp_contact(nop1/numbeads1+nop2/numbeads2,nop1/numbeads1+nop2/numbeads2),dimer(nop1/numbeads1+nop2/numbeads2,nop1/numbeads1+nop2/numbeads2),hp_dimer(nop1/numbeads1+nop2/numbeads2,nop1/numbeads1+nop2/numbeads2)		
        integer pep_nbr(nop1/numbeads1+nop2/numbeads2,10),pep_nbr_num(nop1/numbeads1+nop2/numbeads2),sheet_identity(nop1/numbeads1+nop2/numbeads2),pep_num(nop1/numbeads1+nop2/numbeads2)
		integer sheet_nbr(nop1/numbeads1+nop2/numbeads2,10),sheet_nbr_num(nop1/numbeads1+nop2/numbeads2),fibril_identity(nop1/numbeads1+nop2/numbeads2),num_sheet(nop1/numbeads1+nop2/numbeads2)
        integer sheet_num,pep_count(nop1/numbeads1+nop2/numbeads2),temp_pep_in_sheet(nop1/numbeads1+nop2/numbeads2,nop1/numbeads1+nop2/numbeads2)
        integer fibril_num,sheet_count(nop1/numbeads1+nop2/numbeads2),temp_sheet_in_fibril(nop1/numbeads1+nop2/numbeads2,nop1/numbeads1+nop2/numbeads2)
        integer target_fibril,target_pep_num_in_fibril,move(nop1/numbeads1+nop2/numbeads2),move2(nop1/numbeads1+nop2/numbeads2)
		integer pep_in_sheet(nop1/numbeads1+nop2/numbeads2,nop1/numbeads1+nop2/numbeads2),num_pep_in_sheet(nop1/numbeads1+nop2/numbeads2),sheet_in_fibril(nop1/numbeads1+nop2/numbeads2,nop1/numbeads1+nop2/numbeads2)
		integer pep_in_fibril(nop1/numbeads1+nop2/numbeads2,nop1/numbeads1+nop2/numbeads2),num_pep_in_fibril(nop1/numbeads1+nop2/numbeads2),link_fibril(nop1/numbeads1+nop2/numbeads2,nop1/numbeads1+nop2/numbeads2)   
		integer final_pep_in_sheet(nop1/numbeads1+nop2/numbeads2,nop1/numbeads1+nop2/numbeads2),final_num_pep_in_sheet(nop1/numbeads1+nop2/numbeads2)
		integer final_sheet_in_fibril(nop1/numbeads1+nop2/numbeads2,nop1/numbeads1+nop2/numbeads2),final_num_sheet_in_fibril(nop1/numbeads1+nop2/numbeads2),final_num_pep_in_fibril(nop1/numbeads1+nop2/numbeads2)
		integer final_fibril_num,final_sheet_num,final_sheet_identity(nop1/numbeads1+nop2/numbeads2),final_fibril_identity(nop1/numbeads1+nop2/numbeads2),final_pep_in_fibril(nop1/numbeads1+nop2/numbeads2,nop1/numbeads1+nop2/numbeads2)
        integer force_count,condition1,condition2,fibril_list(nop1/numbeads1+nop2/numbeads2),temp_fibril_list(nop1/numbeads1+nop2/numbeads2)
		logical previous_pep_in_fibril(nop1/numbeads1+nop2/numbeads2)
       end module global
