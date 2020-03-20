# nbodykit parallel issue

mpich from anaconda cannot be used with nbodykit.
Each process is independent and do not communicate.
So we use the system mpi as suggested by Yu.
Only intel mpi and openmpi unfortunately.

Well the computing environment on SDSC popeye sucks.

## email from Yu
Here is one possible way that came into my mind:

Install from bccp, remove mpich forcefully, then pip install pfft-python, mpsort and mpi4py,
hopefully they will use the system mpi implementation.

Or port the nersc directory under conda-channel-bccp there (more work tough).

## mpich to openmpi2
conda activate nbodykit
conda remove --force mpich mpi4py pfft-python mpsort
module load gcc openmpi2
pip install --force-reinstall --no-deps --log pip.log mpi4py
pip install --force-reinstall --no-deps --log pip.log pfft-python mpsort

## openmpi4 complaint
By default, for Open MPI 4.0 and later, infiniband ports on a device
are not used by default.  The intent is to use UCX for these devices.
You can override this policy by setting the btl_openib_allow_ib MCA parameter
to true.

## openmpi3 and openmpi4 glibc version issue
GLIBC_2.27 not found

## nbodykit parallel environment
On SDSC popeye, it requires both conda and module

In slurm script:
source $HOME/anaconda/bin/activate nbodykit
module load gcc openmpi2
