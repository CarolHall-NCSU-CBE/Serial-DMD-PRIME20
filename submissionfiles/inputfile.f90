module allinone
	
	implicit none

	character(len=31) pep1, pep2
	integer nb1,numbeads1,chnln1,nc,nopwg1,nop1
	integer nb2,numbeads2,nb4,chnln2,nc2,nopwg2,nop2
	real*8 boxlength
	real simtemp

	parameter(pep1 = "GVLYVGS")
	parameter(pep2 = "GVLYVGS")
    parameter(nb1=28)          		! number of beads
    parameter(numbeads1=26)         ! number of beads without glycines
    parameter(chnln1=7)       		! chain length
    parameter(nc=24)           		! chain number	
    parameter(nopwg1=nb1*nc)      	! total number of beads
    parameter(nop1=numbeads1*nc)    ! total number of beads minus glycines
    parameter(nb2=28)          		! number of beads
    parameter(numbeads2=26)         ! number of beads without glycines
    parameter(chnln2=7)       		! chain length
    parameter(nc2=24)           	! chain number	
    parameter(nopwg2=nb2*nc2)      	! total number of beads
    parameter(nop2=numbeads2*nc2)   ! total number of beads minus glycines
	parameter(boxlength=159.0D0)	! length of the simulation box in Angstrom
	parameter(simtemp=310.0D0)		! simulation temperature in Kelvin

!	If any value on the second column of output.out is larger than the 0.01 tolerance, adjust
!	d1 d2 d3 values in sidechain adjustment section by increasing the two parameters below
	real*8 dadjust1, dadjust2
	parameter(dadjust1=3.0)
	parameter(dadjust2=3.0)

end module allinone