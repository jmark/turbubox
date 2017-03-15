#!/bin/bash -l

## ---------------------------
## partition gassner-exclusive
##
##                   min max
## Number of nodes   1    4
## Number of sockets 1    2 
## Number of cores   1   12
## Memory per node --  64gb (61 usable)
## Total wall time --  - h
## ---------------------------

#SBATCH --partition=gassner-exclusive
#SBATCH --account=AGWalch
#SBATCH --mail-user=jmarker2@uni-koeln.de
#SBATCH --mail-type=all

#SBATCH --cpus-per-task=1
#SBATCH --nodes=2
#SBATCh --ntasks-per-node=12
#SBATCH --ntasks=24
#SBATCH --mem=40gb
#SBATCH --time=05:00:00

NTASK="$SLURM_NTASKS" "$@"
