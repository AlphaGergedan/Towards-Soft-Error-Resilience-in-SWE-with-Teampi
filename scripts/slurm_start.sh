#!/bin/bash
#SBATCH -J SWE
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=cm2
#SBATCH --partition=cm2_std
#SBATCH --qos=cm2_std
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=end
#SBATCH --mail-user=alexander.hoelzl@tum.de
#SBATCH --export=NONE
#SBATCH --time=03:00:00 

module load slurm_setup
module unload intel-mpi
module load ulfm2/4.0.2u1-gcc8
module load netcdf-hdf5-all/4.7_hdf5-1.10-gcc8-serial  

#Parameters of wrapper Size, numSpares, MTBF, numFails
mpiexec -np $SLURM_NTASKS ./teaMPI_multi_node.sh 4300 0 20 0