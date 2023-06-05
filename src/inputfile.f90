module inputreadin
	implicit none
!	If any value on the second column of output.out is larger than the 0.01 tolerance, adjust
!	d1 d2 d3 values in sidechain adjustment section by increasing the two parameters below
	real*8 dadjust1, dadjust2
	!VN: File number, number setting is arbitary
!VN: File unit 6 is screen or output files
!VN: Files that are stored in running directory:
! Input file
	integer, parameter :: file_input = 1
	integer, parameter :: ideninp = 56
	integer, parameter :: iden2inp = 58
	integer, parameter :: inpinfon1 = 16
	integer, parameter :: inpinfon2 = 18

! Check files
	integer, parameter :: checkid = 60
	integer, parameter :: checkcf1 = 21
	integer, parameter :: checkcf2 = 22
	integer, parameter :: checkcf = 19
	integer, parameter :: checkemblem = 20
	integer, parameter :: checkcfout = 23
	integer, parameter :: checkmass = 24
	integer, parameter :: checksumvel = 25
	integer, parameter :: checkvel = 26
	integer, parameter :: checkhb = 27

! Result parameters:
	integer, parameter :: parahp1 = 9
	integer, parameter :: parahp2 = 10
	integer, parameter :: parafs1 = 11
	integer, parameter :: parafs2 = 12
	integer, parameter :: paraid = 13

! Results:
	integer, parameter :: runcf = 30
	integer, parameter :: runlasvel = 31 
	integer, parameter :: runpartner = 32
	integer, parameter :: rune = 33
	integer, parameter :: runrca = 34
	integer, parameter :: runphipsi = 35
	integer, parameter :: runpdb = 36
	integer, parameter :: precf = 50
	integer, parameter :: prelasvel = 51
	integer, parameter :: prepartner = 52
	integer, parameter :: prerca = 54

! Outputs:
	integer, parameter :: fileout = 55

!VN: Files that are stored in source code package:
!VN: Files that are required for initial configuration
! Generic peptide parameters
	integer, parameter :: genericpepx = 2
	integer, parameter :: genericpepy = 3
	integer, parameter :: genericpepz = 4

! Mass and peptide ID file
	integer, parameter :: massfile = 7
! Bond parameters:
	integer, parameter :: parartoall = 14
	integer, parameter :: parasqz = 15
! Force field for configuration:
	integer, parameter :: ffha55a = 17

!VN: Files opens in inputinfo.f
	integer, parameter :: ep19ha55weak = 37
	integer, parameter :: simwellha55a = 47
	integer, parameter :: parahelix = 38
	integer, parameter :: parabundles = 39
	integer, parameter :: pararn = 40
	integer, parameter :: pararc = 41
	integer, parameter :: parapro = 42
	integer, parameter :: parabdln = 43
	integer, parameter :: simrtoall = 44
	integer, parameter :: simsqz = 45
	integer, parameter :: simmass = 46

!VN: All filenames:
	character*4 fname_digits,fname_digits_pre
	common fname_digits,fname_digits_pre
	character*64 fileinput, fileid1, fileid2, fileinfo_n1, fileinfo_n2
	character*64 fileidck, filecf1ck, filecf2ck, filefinalcheck, filemck, filemassck, filesumvelck, filevelck, filehbck
	character*64 filehp1pr, filehp2pr, filefs1pr, filefs2pr, fileidpr
	character*64 filepepx, filepepy, filepepz, filemasscf
	character*64 filecfrtoall, filecfsqz, filecfha55
	character*64 fileep19w, filesimha55, filehelix, filebundle, filerarn, filerarc, filerapro, filebdln, filesimrall, filesimsqz, fillsimmass


	character(len=31) pep1, pep2
	integer nb1,numbeads1,nopwg1,nop1
	integer nb2,numbeads2,nopwg2,nop2
	!integer nb1p,nb21p,chnln1p
	integer chnln1,nc,numgly1
	integer chnln2,nc2,numgly2
	integer simcoll, numsim
	real*8 boxlength
	real simtemp
	real maxtemp, mintemp
	character*517 ::  path, mydir, rundir
	integer realpath, annealingsteps, simsteps, stepcount, newold, constep
	logical :: back=.true.
	real, allocatable :: temps(:)
	integer, allocatable :: collset(:)

end module inputreadin