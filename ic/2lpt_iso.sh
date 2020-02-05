#!/bin/bash

#PBS -S /bin/bash

#PBS -q mini2
#PBS -l nodes=4:ppn=28
#PBS -l walltime=14:00:00:00
#PBS -j oe
#PBS -o stdout
#PBS -m p
#PBS -M kazuyuki.akitsu@ipmu.jp

#PBS -N ic

cd /work/kazuyuki.akitsu/csit/ic/
mpirun -np 112 ./2LPTic ./masaki_params/masaki_iso.param > ./stdout_iso
