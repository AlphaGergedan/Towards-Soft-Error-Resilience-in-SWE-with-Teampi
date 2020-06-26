#!/bin/bash
#SBATCH -J SWE
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=cm2
#SBATCH --partition=cm2_std
#SBATCH --qos=cm2_std
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=end
#SBATCH --mail-user=alexander.hoelzl@tum.de
#SBATCH --export=NONE
#SBATCH --time=01:30:00 

module load slurm_setup
module unload intel-mpi
module load ulfm2/4.0.2u1-gcc8
module load netcdf-hdf5-all/4.7_hdf5-1.10-gcc8-serial  

export OMPI_MCA_mpi_ft_detector_period=8
export OMPI_MCA_mpi_ft_detector_timeout=24

#Parameters of wrapper Size, numSpares, MTBF, numFails
srun -N $SLURM_NTASKS --ntasks-per-node $SLURM_NTASKS_PER_NODE ./teaMPI_multi_node.sh 4200 0 30 0
srun -N $SLURM_NTASKS --ntasks-per-node $SLURM_NTASKS_PER_NODE ./teaMPI_multi_node.sh 4200 2 30 0 
srun -N $SLURM_NTASKS --ntasks-per-node $SLURM_NTASKS_PER_NODE ./teaMPI_multi_node.sh 4200 4 30 0 
srun -N $SLURM_NTASKS --ntasks-per-node $SLURM_NTASKS_PER_NODE ./teaMPI_multi_node.sh 4200 6 30 0 
srun -N $SLURM_NTASKS --ntasks-per-node $SLURM_NTASKS_PER_NODE ./teaMPI_multi_node.sh 4200 8 30 0 
