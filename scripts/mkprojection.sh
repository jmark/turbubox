#!/usr/bin/env bash

FILES=~/builds/turbulence3D/output/driventurb_3d_hdf5_plt_cnt_0*

for F in $FILES
do
    BASE=$(basename "$F")
    echo $BASE
    ./main "$F" > ./csv/$BASE
done
