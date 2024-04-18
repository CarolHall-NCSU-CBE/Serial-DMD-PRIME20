subroutine readbptnr
	use global
	use inputreadin
	implicit none
	
	integer(kind=k16) collbptnr
	integer :: i,j,argcount, k,l,round, io,annealid,count,ii,jj
	integer izero, start, end, ioerr,collanneal
	integer ichar, len
	real sumanneal,annealtemp,presum
	logical exist_flag
	character(len=2) :: itochar
	character*64 input1
	izero = ichar('0')
	annealid = 0
	sumanneal = 0
	bptnr = 0
	inquire(file = 'inputs/annealtemp_1', exist = exist_flag)
	do annealid = 1,100
		write(itochar,'(i1)') annealid
		itochar = adjustl(itochar)
		input1 = 'inputs/annealtemp_'//itochar
		open(268, file = input1, status = 'old',iostat=ioerr)
		if (ioerr /= 0) then
			go to 103
		endif
		read(268,*) annealtemp
		read(268,*) collanneal
		presum = sumanneal
		sumanneal = sumanneal + real(collanneal)/1000000000
		if (collinbill.le.sumanneal) then
			go to 103
		endif
	enddo
103	close(268)
	if (collinbill .le. sumanneal) then
		i = annealid
	else
		i = 9 + int((collinbill-sumanneal)/(simcoll*1E-9))
	endif
	fname_digits = char(i/1000+izero)//char(mod(i,1000)/100+izero)//char(mod(i,100)/10+izero)//char(mod(i,10)+izero)
	open(runpartner,file = 'results/run'//fname_digits//'.bptnr',status = 'old',form='unformatted')
	do while (.true.)
		read(runpartner) collbptnr,bptnr
		if (abs(collinbill- int((real(collbptnr)/(1E9)+presum)*10000.0)*1E-4) .le. 0.0005) go to 104
		if((collbptnr.gt.0).and.(mod(collbptnr,10000000).eq.0)) goto 104
	enddo
104	close(runpartner)


end subroutine