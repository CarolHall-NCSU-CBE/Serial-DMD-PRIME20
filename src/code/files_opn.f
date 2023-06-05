! This subroutine is to open new files to record results or to read parameters from existing files for simulation step
	subroutine file_opener

#include "def.h"

	use global
	use inputreadin

	implicit none

	integer i, izero
	character*64 filename, zfilename
	logical exist_flag, zexist_flag
	character*1 zero
	character*1 char
	integer ichar, len
	data zero/'0'/
	integer tempwrite,dotpos
	character(len=7) outname
	character*4 simdigit

!VN: For files that are saved in running directory
call chdir(rundir)
! Results:
	izero = ichar('0')
	i = 0
	
	filename = 'results/run'//'0000.energy'
	zfilename = 'results/run'//'0000.energy'//'.gz'
	fname_digits = '0000'
	inquire(file = filename, exist = exist_flag)
	inquire(file = zfilename, exist = zexist_flag)
	do while (exist_flag .or. zexist_flag)		
		i = i + 1
		if (i .gt. 9999) then
			write(0,*) " agf : file name error. too many output files."
			stop
		endif
		fname_digits = char(i/1000+izero)//char(mod(i,1000)/100+izero)//char(mod(i,100)/10+izero)//char(mod(i,10)+izero)
	   	filename = 'results/run'//fname_digits//'.energy'
	   	zfilename = 'results/run'//fname_digits//'.energy'//'.gz'
	   	inquire( file = filename, exist = exist_flag)
	   	inquire( file = zfilename, exist = zexist_flag)
	end do

	i = i -1
	fname_digits_pre = char(i/1000+izero)//char(mod(i,1000)/100+izero)//char(mod(i,100)/10+izero)//char(mod(i,10)+izero)
	open(rune, file=filename, status = 'new')
	filename = 'results/run'//fname_digits//'.config'
	open(unit=runcf,file=filename,status='unknown', form='unformatted',position='append')
	filename = 'results/run'//fname_digits//'.bptnr'
	open(unit=runpartner,file=filename,status='unknown', form='unformatted',position='append')
	filename = 'results/run'//fname_digits//'.lastvel'
	open(unit=runlasvel,file=filename,form='unformatted')
	filename = 'results/run'//fname_digits//'.pdb'
      	open(unit=runpdb, file=filename)

#ifdef write_phipsi
	filename = 'results/run'//fname_digits//'.phipsi'
	open(runphipsi, file=filename, status = 'unknown')
#endif

!VN: Open the data files resulted from previous runs:
	open(unit=precf,file='results/run'//fname_digits_pre//'.config',status='unknown',form='unformatted')
	open(unit=prelasvel,file='results/run'//fname_digits_pre//'.lastvel',status='unknown',form='unformatted')
	open(unit=prepartner,file='results/run'//fname_digits_pre//'.bptnr',status='old',form='unformatted')
	open(unit=prerca,file='results/run'//fname_digits_pre//'.rca',status='unknown')
	open(unit=parafs1,file='parameters/firstside1.data',status='unknown')
	open(unit=parafs2,file='parameters/firstside2.data',status='unknown')
	open(unit=paraid,file='parameters/identity.inp',status='unknown')
	open(unit=parahp1,file='parameters/hp1.inp',status='unknown')
	open(unit=parahp2,file='parameters/hp2.inp',status='unknown')

!VN: Open file for general output:
	write(outname,'(f7.3)') setemp/12
	outname=trim(adjustl(outname))
	dotpos = index(outname,'.')
	outname = outname((dotpos+1):len(outname))
	if (newold == 0) then
	if (stepcount .le. annealingsteps) then 
		open(unit=fileout,file='outputs/out0'//outname,status='unknown')
	else
		simdigit = char((i-annealingsteps+1)/1000+izero)//char(mod((i-annealingsteps+1),1000)/100+izero)//char(mod((i-annealingsteps+1),100)/10+izero)//char(mod((i-annealingsteps+1),10)+izero)
		open(unit=fileout,file='outputs/out0'//outname//'_'//simdigit,status='unknown')
	endif
	else
		simdigit = char((i-annealingsteps+1)/1000+izero)//char(mod((i-annealingsteps+1),1000)/100+izero)//char(mod((i-annealingsteps+1),100)/10+izero)//char(mod((i-annealingsteps+1),10)+izero)
		open(unit=fileout,file='outputs/out0'//outname//'_'//simdigit,status='unknown')
	endif	
	


call chdir(mydir)
! VN: Files are saved in code package:
!     	accessing ideal helix as starting configuration (eight 16mer)
      	open(unit=parabundles,file='parameters/bundles.inp',status='unknown')
	open(unit=pararn,file='parameters/rn.data',status='unknown')
	open(unit=pararc,file='parameters/rc.data',status='unknown')
!     	read in sigma, welldia, and epsilon arrays
      	open(unit=parapro,file='parameters/protein.data',status='unknown')
	open(unit=parabdln,file='parameters/bdln.inp',status='unknown')
	open(unit=ep19ha55weak,file='parametersep/ep19p_ha55a_weakhp.data',status='unknown')
	open(unit=simrtoall,file='parameters/rcarnrco.data',status='unknown')
      	open(unit=simwellha55a,file='parameters/beadwell_ha55a.data',status='unknown')
	open(unit=simsqz,file='parameters/sqz6to10.data',status='old')
	open(unit=simmass,file='parameters/mass.data',status='unknown')








	

	return

	end
