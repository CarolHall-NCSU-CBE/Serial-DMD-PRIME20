	subroutine resultsfile_opener

#include "def.h"

	use global
	use allinone

	implicit none

	integer i, izero
	character*64 filename, zfilename
	logical exist_flag, zexist_flag
	character*1 zero
	character*1 char
	integer ichar, len
	data zero/'0'/
	
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
	open(11111, file=filename, status = 'new')
	filename = 'results/run'//fname_digits//'.config'
	open(unit=800,file=filename,status='unknown', form='unformatted',position='append')
	close(unit=800)
	filename = 'results/run'//fname_digits//'.bptnr'
	open(unit=801,file=filename,status='unknown', form='unformatted',position='append')
	close(unit=801)

#ifdef write_phipsi
	filename = 'results/run'//fname_digits//'.phipsi'
	open(11112, file=filename, status = 'unknown')
#endif

	return

	end
