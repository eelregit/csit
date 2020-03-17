#!/bin/bash

#SBATCH --job-name=res
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kazuyuki.akitsu@ipmu.jp

#SBATCH --partition=general

#SBATCH --exclusive
#SBATCH --nodes=5
#SBATCH --ntasks=240
#SBATCH --mem-per-cpu=16G
#SBATCH --time=3-00:00:00

#SBATCH --chdir=/home/yinli/csit/analysis/

hostname; pwd; date

srun -n $SLURM_NTASKS --mpi-pmi2 python ./response_treepm.py

date
