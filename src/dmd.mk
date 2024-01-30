##MYDMD_PATH = "pwd"

##PROG = DMDPRIME20
## Define Compiler
CC = ifort

## Module file
MODULE = code/inputfile.f90

## Define which code to compile
SOURCE = code/main.F90

## Define array checking flags
OPFLAGS = -Os -xAVX -no-prec-div -r8 -arch host -align dcommons -g -traceback -Dchaptype=1 -Dnumbin=4000 -Dn_wrap=2 -Dcanon -Drunr -Dwrite_phipsi -module code/ -o
## Define executable file
OPEXEC = DMDPRIME20
all:
	$(CC) $(OPFLAGS) $(OPEXEC) $(MODULE) $(SOURCE) 
clean:
	rm -rf code/*.mod
	rm -rf *.o
