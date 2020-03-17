#!/bin/bash

#SBATCH --job-name=ic
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kazuyuki.akitsu@ipmu.jp

#SBATCH --partition=general

#SBATCH --exclusive
#SBATCH --nodes=5
#SBATCH --ntasks=240
#SBATCH --time=2-00:00:00
#SBATCH --exclude=pcn-5-[55-60]
#SBATCH --exclude=pcn-5-[65-69]

#SBATCH --chdir=/home/yinli/csit/ic/

module load gcc openmpi2 lib/fftw2/2.1.5-openmpi2 lib/gsl

hostname; pwd; date

srun ./2LPTic ./planck2015_params/ppm.param

date
