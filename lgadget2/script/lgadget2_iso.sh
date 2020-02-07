#!/bin/bash

#SBATCH --partition=general
#SBATCH --nodes=5
#SBATCH --ntasks=240
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --output=stdout_iso
#SBATCH --workdir=/home/yinli/csit/lgadget2/
#SBATCH --mail-user=kazuyuki.akitsu@ipmu.jp
#SBATCH --mail-type=END
#SBATCH --job-name=nb

/usr/mpi/gcc/openmpi-2.1.2-hfi/bin/mpirun -np 240 ./L-Gadget2 ./masaki_params/params_iso.txt
