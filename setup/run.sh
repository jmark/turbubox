#!/usr/bin/env bash

exec /opt/mpich/bin/mpirun -n $(nproc) ./flash4 | tee stdout+stderr.log
