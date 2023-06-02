##MYDMD_PATH = "pwd"

##PROG = DMDPRIME20
## Define Compiler
CC = ifort

## Module file
MODULE = inputfile.f90

## Define which code to compile
SOURCE = main.F90

## Define array checking flags
CCFLAGS = -Os -xAVX -no-prec-div -r8 -arch host -align dcommons -g -traceback -Dchaptype=1 -Dnumbin=4000 -Dn_wrap=2 -Dcanon -Drunr -Dwrite_phipsi -o

## Define executable file
CCEXEC = ../DMDPRIME20
all:
	$(CC) $(CCFLAGS) $(CCEXEC) $(MODULE) $(SOURCE) 
clean:
	rm *.mod
	rm *.o
