#!/bin/bash -l

# # =========================================================================== #
# # PARTITION MPI
# 
# #SBATCH --partition=mpi
# 
# #SBATCH --cpus-per-task=1
# #SBATCH --nodes=8
# #SBATCh --ntasks-per-node=8
# #SBATCH --ntasks=64
# 
# #SBATCH --mem=24gb

# =========================================================================== #
# PARTITION GASSNER-EXCLUSIVE

#SBATCH --partition=gassner-exclusive

#SBATCH --cpus-per-task=1		# dito
#SBATCH --nodes=3			    # number of nodes needed for the job
#SBATCh --ntasks-per-node=28	# *maximum* number of tasks per node
#SBATCH --ntasks=64			    # number of tasks *actually* needed
#SBATCH --mem=48gb			    # per node memory limit

# With setup above, we have 84 cores at our disposal, but we
# only need 64 tasks due to restrictions of the simulation.
#
# We have 4**3 blocks consisting of 32**3 cells, thus (4 * 32)**3 = 128**3 cells in total.
# These 4**3 = 64 blocks are distributed to all 64 tasks, resp. cores.
# Hence, we waste (84 - 64)/84 = 24% of idle computing power.
#
# Above memory assumption is just an educated guess...
# TODO: Ask hpc support what maximum memory per node is available.

# =========================================================================== #
# GENERAL SETUP

#SBATCH --time=24:00:00			# maximum wall time
# Educated guess...

#SBATCH --account=AGWalch

#SBATCH --job-name=flash.cgs.128.es

#SBATCH --mail-user=jmarker2@uni-koeln.de
#SBATCH --mail-type=all

# =========================================================================== #
# RUNTIME SETUP

module load hdf5
module load intel
module load intelmpi

srun -n $SLURM_NTASKS \
    /usr/bin/time -f "%E real,%U user,%S sys, %D unshared ram" -o runtime.txt \
    ./flash4 flash.par
