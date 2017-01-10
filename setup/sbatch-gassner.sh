#!/bin/bash -l

## ---------------------------
## partition gassner-exclusive
##
##                 min max
## Number of nodes 1   2
## Number of cores 1   16
## Memory per node --  22gb
## Total wall time --  1 h
## ---------------------------

#SBATCH --partition=gassner-exclusive
#SBATCH --account=AGWalch
#SBATCH --mail-user=jmarker2@uni-koeln.de
#SBATCH --mail-type=all

#SBATCH --cpus-per-task=1
#SBATCH --nodes=3
#SBATCh --ntasks-per-node=24
#SBATCH --ntasks=72
#SBATCH --mem=56gb
#SBATCH --time=05:00:00

NTASK="$SLURM_NTASKS" eval "$@"
