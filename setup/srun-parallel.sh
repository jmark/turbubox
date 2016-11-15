#!/bin/sh

set -eu
parallel -N 1 --delay 0.2 -j $SLURM_NTASKS srun -n1 -N1 --exclusive "$@"
