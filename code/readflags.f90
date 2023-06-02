      Logical Function readflags(fileunit,flag)


      IMPLICIT NONE

      Integer :: fileunit
      Character(Len=*) :: flag
      Character(Len=175) :: string

      ! Read string from input file
      Read(fileunit,'(A)') string

      If (flag .eq. string) Then
      	readflags = .true.
      Else
      	readflags = .false.
      Endif

      Return
	End Function



