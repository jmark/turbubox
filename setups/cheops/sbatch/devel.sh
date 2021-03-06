#!/bin/bash --login

## --------------------------
## partion         devel
##
##                 min max
## Number of nodes 1   2
## Number of cores 1   16
## Memory per node --  22gb
## Total wall time --  1 h
## --------------------------

#SBATCH --partition=devel
#SBATCH --account=AGWalch
#SBATCH --mail-user=jmarker2@uni-koeln.de
#SBATCH --mail-type=all

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCh --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=22gb
#SBATCH --time=01:00:00

#SBATCH --parsable

NTASK="$SLURM_NTASKS" "$@"
