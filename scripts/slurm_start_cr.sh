#!/bin/bash
#SBATCH -J SWE
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=end
#SBATCH --mail-user=alexander.hoelzl@tum.de
#SBATCH --export=NONE
#SBATCH --time=02:20:00 

module load slurm_setup
module unload intel-mpi
module load ulfm2/4.0.2u1-gcc8
module load netcdf-hdf5-all/4.7_hdf5-1.10-gcc8-serial  

export OMPI_MCA_mpi_ft_detector_period=20
export OMPI_MCA_mpi_ft_detector_timeout=60
#export OMPI_MCA_mpi_ft_detector_thread=true

export OMP_NUM_THREADS=$((28 / $SLURM_NTASKS_PER_NODE))

NUM_FAILURES=0

for i in $(seq 0 $NUM_FAILURES); do 
	FAIL=$(( i <= $NUM_FAILURES ))
	INIT_RUN=$(( i == 0 )) 
	srun -N $SLURM_JOB_NUM_NODES --ntasks-per-node $SLURM_NTASKS_PER_NODE ./checkpoint_test.sh 3500 30 $FAIL 30 $INIT_RUN
done

