#!/bin/bash

mkdir results
mkdir outputs
mkdir parameters
mkdir inputs
mkdir checks

# Generate initial configuration for new simulation:
../../src_tempaddtoenergyfile/initconfig

#Annealing:	
for i in {1..9}
do ../../src_tempaddtoenergyfile/DMDPRIME20 < inputs/annealtemp_$i > outputs/out_annealtemp_$i
done

#DMD simulations
for i in {1..100}
do ../../src_tempaddtoenergyfile/DMDPRIME20 < inputs/simtemp > outputs/out_simtemp_$i
done