#!/bin/bash

#SBATCH --job-name=nb
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kazuyuki.akitsu@ipmu.jp

#SBATCH --partition=general

#SBATCH --exclusive
#SBATCH --nodes=5
#SBATCH --ntasks=240
#SBATCH --time=2-00:00:00

#SBATCH --chdir=/home/yinli/csit/lgadget2

module load gcc openmpi lib/fftw2/2.1.5-openmpi1 lib/gsl

hostname; pwd; date

srun --mpi=pmi2 ./L-Gadget2 ./planck2015_params/params_iso.txt

date
