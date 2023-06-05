!     subroutine phipsi creates histogram for use in calculating phi-psi angles
!     note, phi,psi bins are fcns of i ... need to change to do 
!     more than one chain

      subroutine phipsi()

#include "def.h"

      use global
	use inputreadin

      implicit none

      integer pos_n,pos_ca,pos_c,pos_n_p1,pos_c_m1,i,j,k,n_phi_shl,n_psi_shl
      real*8 xn,yn,zn,xca,yca,zca,xc,yc,zc,xn_p1,yn_p1,zn_p1,xc_m1,yc_m1,zc_m1
      real*8 x_nca,y_nca,z_nca,dnca,x1_pri,y1_pri,z1_pri,x_cac,y_cac,z_cac,dcac
      real*8 cosangle,angle,dcaq,xq,yq,zq,x_qc,y_qc,z_qc,dcq,xr,yr,zr,dcca
      real*8 x2_pri,y2_pri,z2_pri,x3_pri,y3_pri,z3_pri,xdiff,ydiff,zdiff
      real*8 xcm1_new,ycm1_new,zcm1_new,xs,ys,zs,x_scm1,y_scm1,z_scm1,dscm1
      real*8 x_can,y_can,z_can,dcan,x_snp,y_snp1,z_snp1,dsnp1,x_qn,y_qn,z_qn,dnq
      real*8 x_cca,y_cca,z_cca,xnp1_newz,ynp1_newz,znp1_new,degs,xnp1_new,ynp1_new,x_snp1

      do k=1,nop1/numbeads1

!      cycle through residues 2 through chnln-1 for each chain
       do i=(k-1)*numbeads1+2,(k-1)*numbeads1+chnln1-1
	  
         pos_n=chnln1+i
         xn=sv(1,pos_n)+sv(4,pos_n)*tfalse
         yn=sv(2,pos_n)+sv(5,pos_n)*tfalse
         zn=sv(3,pos_n)+sv(6,pos_n)*tfalse
         
         pos_ca=i
         xca=sv(1,pos_ca)+sv(4,pos_ca)*tfalse
         yca=sv(2,pos_ca)+sv(5,pos_ca)*tfalse
         zca=sv(3,pos_ca)+sv(6,pos_ca)*tfalse
         
         pos_c=2*chnln1+i
         xc=sv(1,pos_c)+sv(4,pos_c)*tfalse
         yc=sv(2,pos_c)+sv(5,pos_c)*tfalse
         zc=sv(3,pos_c)+sv(6,pos_c)*tfalse
         
         pos_n_p1=chnln1+i+1
         xn_p1=sv(1,pos_n_p1)+sv(4,pos_n_p1)*tfalse
         yn_p1=sv(2,pos_n_p1)+sv(5,pos_n_p1)*tfalse
         zn_p1=sv(3,pos_n_p1)+sv(6,pos_n_p1)*tfalse
         
         pos_c_m1=2*chnln1+i-1
         xc_m1=sv(1,pos_c_m1)+sv(4,pos_c_m1)*tfalse
         yc_m1=sv(2,pos_c_m1)+sv(5,pos_c_m1)*tfalse
         zc_m1=sv(3,pos_c_m1)+sv(6,pos_c_m1)*tfalse

!        phi

!        find unit vector from ni to cai 
!        this is x-axis in new frame of reference
         x_nca=xca-xn
         y_nca=yca-yn
         z_nca=zca-zn
	 x_nca=x_nca-dnint(x_nca)
	 y_nca=y_nca-dnint(y_nca)
	 z_nca=z_nca-dnint(z_nca)
         dnca=dsqrt(x_nca*x_nca+y_nca*y_nca+z_nca*z_nca)
         x1_pri=x_nca/dnca
         y1_pri=y_nca/dnca
         z1_pri=z_nca/dnca

!        find unit vector from cai to ci
         x_cac=xc-xca
         y_cac=yc-yca
         z_cac=zc-zca
	 x_cac=x_cac-dnint(x_cac)
	 y_cac=y_cac-dnint(y_cac)
	 z_cac=z_cac-dnint(z_cac)
         dcac=dsqrt(x_cac*x_cac+y_cac*y_cac+z_cac*z_cac)
         x_cac=x_cac/dcac
         y_cac=y_cac/dcac
         z_cac=z_cac/dcac

!        find angle between the two vectors above
         cosangle=x1_pri*x_cac+y1_pri*y_cac+z1_pri*z_cac

         if (cosangle.gt.1.d0) cosangle=1.d0
         if (cosangle.lt.-1.d0) cosangle=-1.d0

         angle=acos(cosangle)
!        degs=180*angle/pi
               
!        find distance btwn cai and q
!        q is the point along line made by ni and cia where line ci-q
!        intersects ni-cia at a right angle
         dcaq=cosangle*dcac

!        find coordinates of q
         xq=xca+x1_pri*dcaq
         yq=yca+y1_pri*dcaq
         zq=zca+z1_pri*dcaq

!        find vector from q to ci
         x_qc=xc-xq
         y_qc=yc-yq
         z_qc=zc-zq
	 x_qc=x_qc-dnint(x_qc)
	 y_qc=y_qc-dnint(y_qc)
	 z_qc=z_qc-dnint(z_qc)
         dcq=dsqrt(dcac*dcac-dcaq*dcaq)
         x_qc=x_qc/dcq
         y_qc=y_qc/dcq
         z_qc=z_qc/dcq

!        find r (the end of the vector made by sliding q - ci onto ni
!        so that q lies on ni and ci lies on r)
         xr=xc-x1_pri*(dnca+dcaq)
         yr=yc-y1_pri*(dnca+dcaq)
         zr=zc-z1_pri*(dnca+dcaq)

!        find unit vector from ni to r
!        this is y-axis in new frame of reference
	 x2_pri=xr-xn
         y2_pri=yr-yn
         z2_pri=zr-zn
	 x2_pri=x2_pri-dnint(x2_pri)
	 y2_pri=y2_pri-dnint(y2_pri)
	 z2_pri=z2_pri-dnint(z2_pri)
         x2_pri=x2_pri/dcq
         y2_pri=y2_pri/dcq
         z2_pri=z2_pri/dcq

!        find third unit vector by taking rh-cross product
         x3_pri=y1_pri*z2_pri-z1_pri*y2_pri
         y3_pri=z1_pri*x2_pri-x1_pri*z2_pri
         z3_pri=x1_pri*y2_pri-y1_pri*x2_pri

!        find coords of c i-1 in new frame of reference
         xdiff=xc_m1-xn
         ydiff=yc_m1-yn
         zdiff=zc_m1-zn
	 xdiff=xdiff-dnint(xdiff)
	 ydiff=ydiff-dnint(ydiff)
	 zdiff=zdiff-dnint(zdiff)
         xcm1_new=x1_pri*xdiff+y1_pri*ydiff+z1_pri*zdiff
         ycm1_new=x2_pri*xdiff+y2_pri*ydiff+z2_pri*zdiff
         zcm1_new=x3_pri*xdiff+y3_pri*ydiff+z3_pri*zdiff

!        define s 
!        sis the point along line made by ni and cia where line ciminus1-s
!        intersects ni-cia at a right angle
         xs=xcm1_new
         ys=0.d0
         zs=0.d0

!        find the angle s-ciminus1 makes with nr (both in new frm of ref)
!        in new frm of ref, nr is (0,1,0)
         x_scm1=0.d0
         y_scm1=ycm1_new
         z_scm1=zcm1_new
         dscm1=dsqrt(ycm1_new*ycm1_new+zcm1_new*zcm1_new)
         y_scm1=ycm1_new/dscm1
         z_scm1=zcm1_new/dscm1

         cosangle=y_scm1
         if (cosangle.gt.1.d0) cosangle=1.d0
         if (cosangle.lt.-1.d0) cosangle=-1.d0
         angle=acos(cosangle)
         degs=180.d0*angle/pi

         if (z_scm1.gt.0.d0) degs=degs*-1.d0

         if (degs.lt.-180.d0) degs=-180.d0
         if (degs.gt.180.d0) degs=180.d0

         degs=degs+180.d0

         n_phi_shl=dint(degs)
        
!        psi

!        find unit vector from ci to cai 
!        this is x-axis in new frame of reference
         x_cca=xca-xc
         y_cca=yca-yc
         z_cca=zca-zc
	 x_cca=x_cca-dnint(x_cca)
	 y_cca=y_cca-dnint(y_cca)
	 z_cca=z_cca-dnint(z_cca)
         dcca=dsqrt(x_cca*x_cca+y_cca*y_cca+z_cca*z_cca)
         x1_pri=x_cca/dcca
         y1_pri=y_cca/dcca
         z1_pri=z_cca/dcca

!        find unit vector from cai to ni
         x_can=xn-xca
         y_can=yn-yca
         z_can=zn-zca
	 x_can=x_can-dnint(x_can)
	 y_can=y_can-dnint(y_can)
	 z_can=z_can-dnint(z_can)
         dcan=dsqrt(x_can*x_can+y_can*y_can+z_can*z_can)
         x_can=x_can/dcan
         y_can=y_can/dcan
         z_can=z_can/dcan

!        find angle between the two vectors above
         cosangle=x1_pri*x_can+y1_pri*y_can+z1_pri*z_can
         if (cosangle.gt.1.d0) cosangle=1.d0
         if (cosangle.lt.-1.d0) cosangle=-1.d0
         angle=acos(cosangle)
!        degs=180*angle/pi
               
!        find distance btwn cai and q
!        q is the point along line made by ci and cia where line ni-q
!        intersects ci-cia at a right angle
         dcaq=cosangle*dcan

!        find coordinates of q
         xq=xca+x1_pri*dcaq
         yq=yca+y1_pri*dcaq
         zq=zca+z1_pri*dcaq

!        find vector from q to ni
         x_qn=xn-xq
         y_qn=yn-yq
         z_qn=zn-zq
	 x_qn=x_qn-dnint(x_qn)
	 y_qn=y_qn-dnint(y_qn)
	 z_qn=z_qn-dnint(z_qn)
         dnq=dsqrt(dcan*dcan-dcaq*dcaq)
         x_qn=x_qn/dnq
         y_qn=y_qn/dnq
         z_qn=z_qn/dnq

!        find r (the end of the vector made by sliding q - ni onto ci
!        so that q lies on ci and ni lies on r)
         xr=xn-x1_pri*(dcca+dcaq)
         yr=yn-y1_pri*(dcca+dcaq)
         zr=zn-z1_pri*(dcca+dcaq)

!        find unit vector from ci to r this is y-axis in new frame of reference
         x2_pri=xr-xc
         y2_pri=yr-yc
         z2_pri=zr-zc
	 x2_pri=x2_pri-dnint(x2_pri)
	 y2_pri=y2_pri-dnint(y2_pri)
	 z2_pri=z2_pri-dnint(z2_pri)
         x2_pri=x2_pri/dnq
         y2_pri=y2_pri/dnq
         z2_pri=z2_pri/dnq

!        find third unit vector by taking rh-cross product
         x3_pri=y1_pri*z2_pri-z1_pri*y2_pri
         y3_pri=z1_pri*x2_pri-x1_pri*z2_pri
         z3_pri=x1_pri*y2_pri-y1_pri*x2_pri

!        find coords of n i+1 in new frame of reference
         xdiff=xn_p1-xc
         ydiff=yn_p1-yc
         zdiff=zn_p1-zc
	 xdiff=xdiff-dnint(xdiff)
	 ydiff=ydiff-dnint(ydiff)
	 zdiff=zdiff-dnint(zdiff)
         xnp1_new=x1_pri*xdiff+y1_pri*ydiff+z1_pri*zdiff
         ynp1_new=x2_pri*xdiff+y2_pri*ydiff+z2_pri*zdiff
         znp1_new=x3_pri*xdiff+y3_pri*ydiff+z3_pri*zdiff

!        define s 
!        s is the point along line made by ci and cia where line niplus1-s
!        intersects ci-cia at a right angle
         xs=xnp1_new
         ys=0.d0
         zs=0.d0

!        find the angle s-niplus1 makes with cr (both in new frm of ref)
!        in new frm of ref, cr is (0,1,0)
         x_snp1=0.d0
         y_snp1=ynp1_new
         z_snp1=znp1_new
         dsnp1=dsqrt(ynp1_new*ynp1_new+znp1_new*znp1_new)
         y_snp1=ynp1_new/dsnp1
         z_snp1=znp1_new/dsnp1

         cosangle=y_snp1
         if (cosangle.gt.1.d0) cosangle=1.d0
         if (cosangle.lt.-1.d0) cosangle=-1.d0
         angle=acos(cosangle)
         degs=180.d0*angle/pi

         if (z_snp1.gt.0.d0) degs=degs*-1.d0

         if (degs.lt.-180.d0) degs=-180.d0
         if (degs.gt.180.d0) degs=180.d0

         degs=degs+180.d0
         n_psi_shl=dint(degs)

         res(n_phi_shl,n_psi_shl)=res(n_phi_shl,n_psi_shl)+1

       enddo

      enddo

do k=(nop1/numbeads1)+1,numchains

!      cycle through residues 2 through chnln-1 for each chain
       do i=(k-1)*numbeads2+2,(k-1)*numbeads2+chnln2-1
	  
         pos_n=chnln2+i
         xn=sv(1,pos_n)+sv(4,pos_n)*tfalse
         yn=sv(2,pos_n)+sv(5,pos_n)*tfalse
         zn=sv(3,pos_n)+sv(6,pos_n)*tfalse
         
         pos_ca=i
         xca=sv(1,pos_ca)+sv(4,pos_ca)*tfalse
         yca=sv(2,pos_ca)+sv(5,pos_ca)*tfalse
         zca=sv(3,pos_ca)+sv(6,pos_ca)*tfalse
         
         pos_c=2*chnln2+i
         xc=sv(1,pos_c)+sv(4,pos_c)*tfalse
         yc=sv(2,pos_c)+sv(5,pos_c)*tfalse
         zc=sv(3,pos_c)+sv(6,pos_c)*tfalse
         
         pos_n_p1=chnln2+i+1
         xn_p1=sv(1,pos_n_p1)+sv(4,pos_n_p1)*tfalse
         yn_p1=sv(2,pos_n_p1)+sv(5,pos_n_p1)*tfalse
         zn_p1=sv(3,pos_n_p1)+sv(6,pos_n_p1)*tfalse
         
         pos_c_m1=2*chnln2+i-1
         xc_m1=sv(1,pos_c_m1)+sv(4,pos_c_m1)*tfalse
         yc_m1=sv(2,pos_c_m1)+sv(5,pos_c_m1)*tfalse
         zc_m1=sv(3,pos_c_m1)+sv(6,pos_c_m1)*tfalse

!        phi

!        find unit vector from ni to cai 
!        this is x-axis in new frame of reference
         x_nca=xca-xn
         y_nca=yca-yn
         z_nca=zca-zn
	 x_nca=x_nca-dnint(x_nca)
	 y_nca=y_nca-dnint(y_nca)
	 z_nca=z_nca-dnint(z_nca)
         dnca=dsqrt(x_nca*x_nca+y_nca*y_nca+z_nca*z_nca)
         x1_pri=x_nca/dnca
         y1_pri=y_nca/dnca
         z1_pri=z_nca/dnca

!        find unit vector from cai to ci
         x_cac=xc-xca
         y_cac=yc-yca
         z_cac=zc-zca
	 x_cac=x_cac-dnint(x_cac)
	 y_cac=y_cac-dnint(y_cac)
	 z_cac=z_cac-dnint(z_cac)
         dcac=dsqrt(x_cac*x_cac+y_cac*y_cac+z_cac*z_cac)
         x_cac=x_cac/dcac
         y_cac=y_cac/dcac
         z_cac=z_cac/dcac

!        find angle between the two vectors above
         cosangle=x1_pri*x_cac+y1_pri*y_cac+z1_pri*z_cac

         if (cosangle.gt.1.d0) cosangle=1.d0
         if (cosangle.lt.-1.d0) cosangle=-1.d0

         angle=acos(cosangle)
!        degs=180*angle/pi
               
!        find distance btwn cai and q
!        q is the point along line made by ni and cia where line ci-q
!        intersects ni-cia at a right angle
         dcaq=cosangle*dcac

!        find coordinates of q
         xq=xca+x1_pri*dcaq
         yq=yca+y1_pri*dcaq
         zq=zca+z1_pri*dcaq

!        find vector from q to ci
         x_qc=xc-xq
         y_qc=yc-yq
         z_qc=zc-zq
	 x_qc=x_qc-dnint(x_qc)
	 y_qc=y_qc-dnint(y_qc)
	 z_qc=z_qc-dnint(z_qc)
         dcq=dsqrt(dcac*dcac-dcaq*dcaq)
         x_qc=x_qc/dcq
         y_qc=y_qc/dcq
         z_qc=z_qc/dcq

!        find r (the end of the vector made by sliding q - ci onto ni
!        so that q lies on ni and ci lies on r)
         xr=xc-x1_pri*(dnca+dcaq)
         yr=yc-y1_pri*(dnca+dcaq)
         zr=zc-z1_pri*(dnca+dcaq)

!        find unit vector from ni to r
!        this is y-axis in new frame of reference
	 x2_pri=xr-xn
         y2_pri=yr-yn
         z2_pri=zr-zn
	 x2_pri=x2_pri-dnint(x2_pri)
	 y2_pri=y2_pri-dnint(y2_pri)
	 z2_pri=z2_pri-dnint(z2_pri)
         x2_pri=x2_pri/dcq
         y2_pri=y2_pri/dcq
         z2_pri=z2_pri/dcq

!        find third unit vector by taking rh-cross product
         x3_pri=y1_pri*z2_pri-z1_pri*y2_pri
         y3_pri=z1_pri*x2_pri-x1_pri*z2_pri
         z3_pri=x1_pri*y2_pri-y1_pri*x2_pri

!        find coords of c i-1 in new frame of reference
         xdiff=xc_m1-xn
         ydiff=yc_m1-yn
         zdiff=zc_m1-zn
	 xdiff=xdiff-dnint(xdiff)
	 ydiff=ydiff-dnint(ydiff)
	 zdiff=zdiff-dnint(zdiff)
         xcm1_new=x1_pri*xdiff+y1_pri*ydiff+z1_pri*zdiff
         ycm1_new=x2_pri*xdiff+y2_pri*ydiff+z2_pri*zdiff
         zcm1_new=x3_pri*xdiff+y3_pri*ydiff+z3_pri*zdiff

!        define s 
!        sis the point along line made by ni and cia where line ciminus1-s
!        intersects ni-cia at a right angle
         xs=xcm1_new
         ys=0.d0
         zs=0.d0

!        find the angle s-ciminus1 makes with nr (both in new frm of ref)
!        in new frm of ref, nr is (0,1,0)
         x_scm1=0.d0
         y_scm1=ycm1_new
         z_scm1=zcm1_new
         dscm1=dsqrt(ycm1_new*ycm1_new+zcm1_new*zcm1_new)
         y_scm1=ycm1_new/dscm1
         z_scm1=zcm1_new/dscm1

         cosangle=y_scm1
         if (cosangle.gt.1.d0) cosangle=1.d0
         if (cosangle.lt.-1.d0) cosangle=-1.d0
         angle=acos(cosangle)
         degs=180.d0*angle/pi

         if (z_scm1.gt.0.d0) degs=degs*-1.d0

         if (degs.lt.-180.d0) degs=-180.d0
         if (degs.gt.180.d0) degs=180.d0

         degs=degs+180.d0

         n_phi_shl=dint(degs)
        
!        psi

!        find unit vector from ci to cai 
!        this is x-axis in new frame of reference
         x_cca=xca-xc
         y_cca=yca-yc
         z_cca=zca-zc
	 x_cca=x_cca-dnint(x_cca)
	 y_cca=y_cca-dnint(y_cca)
	 z_cca=z_cca-dnint(z_cca)
         dcca=dsqrt(x_cca*x_cca+y_cca*y_cca+z_cca*z_cca)
         x1_pri=x_cca/dcca
         y1_pri=y_cca/dcca
         z1_pri=z_cca/dcca

!        find unit vector from cai to ni
         x_can=xn-xca
         y_can=yn-yca
         z_can=zn-zca
	 x_can=x_can-dnint(x_can)
	 y_can=y_can-dnint(y_can)
	 z_can=z_can-dnint(z_can)
         dcan=dsqrt(x_can*x_can+y_can*y_can+z_can*z_can)
         x_can=x_can/dcan
         y_can=y_can/dcan
         z_can=z_can/dcan

!        find angle between the two vectors above
         cosangle=x1_pri*x_can+y1_pri*y_can+z1_pri*z_can
         if (cosangle.gt.1.d0) cosangle=1.d0
         if (cosangle.lt.-1.d0) cosangle=-1.d0
         angle=acos(cosangle)
!        degs=180*angle/pi
               
!        find distance btwn cai and q
!        q is the point along line made by ci and cia where line ni-q
!        intersects ci-cia at a right angle
         dcaq=cosangle*dcan

!        find coordinates of q
         xq=xca+x1_pri*dcaq
         yq=yca+y1_pri*dcaq
         zq=zca+z1_pri*dcaq

!        find vector from q to ni
         x_qn=xn-xq
         y_qn=yn-yq
         z_qn=zn-zq
	 x_qn=x_qn-dnint(x_qn)
	 y_qn=y_qn-dnint(y_qn)
	 z_qn=z_qn-dnint(z_qn)
         dnq=dsqrt(dcan*dcan-dcaq*dcaq)
         x_qn=x_qn/dnq
         y_qn=y_qn/dnq
         z_qn=z_qn/dnq

!        find r (the end of the vector made by sliding q - ni onto ci
!        so that q lies on ci and ni lies on r)
         xr=xn-x1_pri*(dcca+dcaq)
         yr=yn-y1_pri*(dcca+dcaq)
         zr=zn-z1_pri*(dcca+dcaq)

!        find unit vector from ci to r this is y-axis in new frame of reference
         x2_pri=xr-xc
         y2_pri=yr-yc
         z2_pri=zr-zc
	 x2_pri=x2_pri-dnint(x2_pri)
	 y2_pri=y2_pri-dnint(y2_pri)
	 z2_pri=z2_pri-dnint(z2_pri)
         x2_pri=x2_pri/dnq
         y2_pri=y2_pri/dnq
         z2_pri=z2_pri/dnq

!        find third unit vector by taking rh-cross product
         x3_pri=y1_pri*z2_pri-z1_pri*y2_pri
         y3_pri=z1_pri*x2_pri-x1_pri*z2_pri
         z3_pri=x1_pri*y2_pri-y1_pri*x2_pri

!        find coords of n i+1 in new frame of reference
         xdiff=xn_p1-xc
         ydiff=yn_p1-yc
         zdiff=zn_p1-zc
	 xdiff=xdiff-dnint(xdiff)
	 ydiff=ydiff-dnint(ydiff)
	 zdiff=zdiff-dnint(zdiff)
         xnp1_new=x1_pri*xdiff+y1_pri*ydiff+z1_pri*zdiff
         ynp1_new=x2_pri*xdiff+y2_pri*ydiff+z2_pri*zdiff
         znp1_new=x3_pri*xdiff+y3_pri*ydiff+z3_pri*zdiff

!        define s 
!        s is the point along line made by ci and cia where line niplus1-s
!        intersects ci-cia at a right angle
         xs=xnp1_new
         ys=0.d0
         zs=0.d0

!        find the angle s-niplus1 makes with cr (both in new frm of ref)
!        in new frm of ref, cr is (0,1,0)
         x_snp1=0.d0
         y_snp1=ynp1_new
         z_snp1=znp1_new
         dsnp1=dsqrt(ynp1_new*ynp1_new+znp1_new*znp1_new)
         y_snp1=ynp1_new/dsnp1
         z_snp1=znp1_new/dsnp1

         cosangle=y_snp1
         if (cosangle.gt.1.d0) cosangle=1.d0
         if (cosangle.lt.-1.d0) cosangle=-1.d0
         angle=acos(cosangle)
         degs=180.d0*angle/pi

         if (z_snp1.gt.0.d0) degs=degs*-1.d0

         if (degs.lt.-180.d0) degs=-180.d0
         if (degs.gt.180.d0) degs=180.d0

         degs=degs+180.d0
         n_psi_shl=dint(degs)

         res(n_phi_shl,n_psi_shl)=res(n_phi_shl,n_psi_shl)+1

       enddo

      enddo


      return

      end

       
