#!/bin/bash

#SBATCH --job-name=fof
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

echo "### conda environment ###"
source $HOME/anaconda/bin/activate nbodykit
conda info

echo "### module environment ###"
module purge
module load slurm
#module load gcc openmpi
module load gcc openmpi2
#module load modules-nix nix/gcc/8.3.0 nix/openmpi3
#module load modules-nix nix/openmpi3
#module load modules-nix nix/gcc nix/openmpi4
module list

srun python3 fof.py

date
