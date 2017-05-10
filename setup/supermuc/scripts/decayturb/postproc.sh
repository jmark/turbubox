#!/bin/bash
~/scripts/parallel/llrun.pl -f ../hw.micro.conf :: :: NPROCS=3 ../plot.sh sim_State_*.h5 | llsubmit -
~/scripts/parallel/llrun.pl -f ../hw.micro.conf :: :: NPROCS=3 ../fft.sh sim_State_*.h5 | llsubmit -
