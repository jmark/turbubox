# Synopsis
    In this folder there are several scripts for wrangling meaningful data out
    of the checkpoint/plots files produced by FLASH.

## py.sh
    Preloader for all python scripts in this directory. Run like so:
        
    $ py.sh ${myscript}.py arg0 arg1 arg2 ...

## evolution.py
    Loop over given files via stdin and condense data to global quantities of
    the simulation, like total kinetic energy, etc...

    $ find /scratch/jmarker2/stirturb/cgs.128.es/amr -name 'flash_hdf5_chk_0*' -type f \
        | sort | py.sh evolution.py | tee result.dat

## powerspectrum.*.py
    Compute powersectrum of given snapshot.
   
    $ py.sh powerspectrum.cython.py \
        /scratch/jmarker2/stirturb/cgs.128.b5/amr/flash_hdf5_chk_0035 \
        | tee result.dat
