#!/bin/bash

#SBATCH --job-name=res
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

srun -n 240 python ./response_pm512.py

date
