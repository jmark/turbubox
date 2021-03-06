#!/bin/bash -l

## --------------------------
## partition       mpi
##
##                 min max
## Number of nodes 2   128
## Number of cores 9   1536
## Memory per node --  46gb
## Total wall time --  360h
## --------------------------

#SBATCH --partition=mpi
#SBATCH --account=AGWalch
#SBATCH --mail-user=jmarker2@uni-koeln.de
#SBATCH --mail-type=all

#SBATCH --cpus-per-task=1
#SBATCH --nodes=8
#SBATCh --ntasks-per-node=8
#SBATCH --ntasks=64
#SBATCH --mem=22gb
#SBATCH --time=24:00:00

module purge
module load hdf5
module load intel
module load intelmpi

srun -n "$SLURM_NTASKS" "$@"
