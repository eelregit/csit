#!/bin/bash

#SBATCH --partition=general
#SBATCH --nodes=6
#SBATCH --ntasks=288
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --output=stdout_ppm
#SBATCH --workdir=/home/yinli/csit/ic/
#SBATCH --mail-user=kazuyuki.akitsu@ipmu.jp
#SBATCH --mail-type=END
#SBATCH --job-name=ic

/usr/mpi/gcc/openmpi-2.1.2-hfi/bin/mpirun -np 288 ./2LPTic ./masaki_params/masaki_ppm.param
