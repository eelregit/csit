#!/bin/bash

#SBATCH --job-name=pow
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kazuyuki.akitsu@ipmu.jp

#SBATCH --partition=general

#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --time=3-00:00:00

#SBATCH --chdir=/home/yinli/csit/analysis/

hostname; pwd; date

source $HOME/anaconda/bin/activate nbodykit

srun python ./power.py

date
