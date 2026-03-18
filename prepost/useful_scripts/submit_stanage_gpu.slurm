#!/bin/bash

#SBATCH --time=00:20:00
##SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --partition=gpu
##SBATCH --qos=short
#SBATCH --qos=gpu
#SBATCH --gres=gpu:4

module use /users/me1blx/SOFTWARE/compilers/nvhpc/modulefiles
module load FFTW.MPI/3.3.10-gompi-2022a
module load nvhpc/22.11
export PATH=$HOME/SOFTWARE/nsightSystems/target-linux-x64:$PATH
export NVHPC=$HOME/SOFTWARE/compilers/nvhpc/Linux_x86_64/22.11
export PATH=$NVHPC/comm_libs/openmpi/openmpi-3.1.5/bin:$PATH
export LD_LIBRARY_PATH=$NVHPC/comm_libs/openmpi/openmpi-3.1.5/lib:$LD_LIBRARY_PATH

mpirun -n 4 bind.sh ../../bin/CHAPSim
