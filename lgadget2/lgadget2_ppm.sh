#!/bin/bash

#PBS -S /bin/bash

#PBS -q mini
#PBS -l nodes=1:ppn=28
#PBS -l walltime=14:00:00:00
#PBS -j oe
#PBS -o stdout
#PBS -m p
#PBS -M kazuyuki.akitsu@ipmu.jp

#PBS -N lgadget2

cd /work/kazuyuki.akitsu/csit/lgadget2/
mpirun -np 28 ./L-Gadget2 ./masaki_params/params_ppm.txt > ./stdout_ppm
