#!/bin/bash

#SBATCH --job-name=fof
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kazuyuki.akitsu@ipmu.jp

#SBATCH --partition=general

#SBATCH --exclusive
#SBATCH --nodes=5
#SBATCH --ntasks=240
#SBATCH --time=3-00:00:00

#SBATCH --chdir=/home/yinli/csit/analysis/

hostname; pwd; date

srun -n $SLURM_NTASKS python ./fof.py

date
