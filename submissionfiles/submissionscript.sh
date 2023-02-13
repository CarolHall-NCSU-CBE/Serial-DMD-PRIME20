#!/bin/bash

# Delete all files of previous simulations if the package is copied over from an old simulation.
cd ../genconfig/
rm -f chninfo*.data
rm -f results/*
rm -f parameters/firstside*.*
rm -f parameters/hp*.*
rm -f parameters/identity.inp
rm -f checks/*
rm -f inputs/id*.*
rm -f seq*.*
cd ../parameters/
rm -f firstside*.*
rm -f hp*.*
rm -f identity.inp
cd ../results
rm -f run*.*
rm -f ../dmd
cd ../code/
rm -f allinone.mod
cd ../temperatures/
rm -f simtemp
cd ../outputs/
rm -f out*

# Generate initial configuration
cd ../submissionfiles/
rm -f allinone.mod
ifort -O3 -o ../genconfig/genconfig inputfile.f90 ../genconfig/scriptinput_gen-SQZ.f90
cd ../genconfig
./genconfig

if [ "$(ls -A results/)" ]
then
	#Copy initial configuration to simulation directory
	cp ../genconfig/results/*.* ../results/
	cp ../genconfig/parameters/firstside*.* ../parameters/
	cp ../genconfig/parameters/hp*.* ../parameters/
	cp ../genconfig/parameters/identity.inp ../parameters/
else
   exit 1 # exit with unspecified error code
fi 

cd ../code

ifort -Os -xAVX -no-prec-div -r8 -arch host -align dcommons -g -traceback -Dchaptype=1 -Dnumbin=4000 -Dn_wrap=2 -Dcanon -Drunr -Dwrite_phipsi -o ../dmd ../submissionfiles/inputfile.f90 main.F90

cd ../
temps='050 045 040 035 030 028 026 024 022' 
for i in $temps
do ./dmd < temperatures/temp_$i > outputs/out_$i wait
done

for i in {1..50}
do ./dmd < temperatures/simtemp > outputs/out_simtemp_$i wait
done





 
