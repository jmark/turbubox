# Synopsis
    In this folder there are several scripts for wrangling meaningful data out
    of the checkpoint/plots files produced by FLASH.

- python/
    Python scripts. These are never run directly. See remarks for bin/ directory.

- bin/
    Preloader for all python scripts in python/ directory. Run like so:
        
    $ ./bin/<script name> arg0 arg1 arg2 ...

    - evolution
        Loop over given files via stdin and condense data to global quantities of
        the simulation, like total kinetic energy, etc...

        $ find <path with checkpoint files> -name 'flash_hdf5_chk_0*' -type f \
            | sort | ./bin/evolution | tee result.dat

    - powerspectrum.*.py
        Compute powersectrum of given snapshot.
       
        $ ./bin/powerspectrum.cython \
            <path with checkpoint files>/flash_hdf5_chk_0XXX | tee result.dat
