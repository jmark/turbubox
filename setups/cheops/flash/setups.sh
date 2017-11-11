#!/bin/sh

set -eu -o pipefail
test -z "$*" && echo "ERROR in $0: No arguments given!" 1>&2 && exit 1

# SOLVER='--with-unit=physics/Hydro/HydroMain/split/PPM'
# SOLVER='--with-unit=physics/Hydro/HydroMain/split/MHD_8Wave'
# SOLVER='--with-unit=physics/Hydro/HydroMain/split/Bouchut3'
# SOLVER='--with-unit=physics/Hydro/HydroMain/split/Bouchut5'
# SOLVER='--with-unit=physics/Hydro/HydroMain/split/ES'

for arg in "$@"; do case "$arg" in

    "KHI 8wave")
        BLOCKSIZE=96

        MKFDIR="$HOME/turbubox/setups/cheops/flash/makefiles/intel"
        RUNDIR="/scratch/jmarker2/flash/KHI/8wave"
        OBJDIR="${RUNDIR}/build"

        mkdir -p "$RUNDIR/chkpt" "$OBJDIR"

        ./setup 'Markert-KHI' -2d +ug -auto -opt -portable \
            '--with-unit=physics/Hydro/HydroMain/split/MHD_8Wave' \
            -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -maxblocks=-1 \
            -site="$MKFDIR" -objdir="$OBJDIR"
    ;;

    "KHI PPM")
        BLOCKSIZE=96

        MKFDIR="$HOME/turbubox/setups/cheops/flash/makefiles/intel"
        RUNDIR="/scratch/jmarker2/flash/KHI/PPM"
        OBJDIR="${RUNDIR}/build"

        mkdir -p "$RUNDIR/chkpt" "$OBJDIR"

        ./setup 'Markert-KHI' -2d +ug -auto -opt -portable \
            '--with-unit=physics/Hydro/HydroMain/split/PPM' \
            -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -maxblocks=-1 \
            -site="$MKFDIR" -objdir="$OBJDIR"
    ;;

    "KHI ES")
        BLOCKSIZE=96
        MKFDIR=$HOME/turbubox/setups/cheops/flash/makefiles/intel
        OBJDIR=/scratch/jmarker2/flash/KHI/ES/build

        mkdir -p "$OBJDIR"

        ./setup 'Markert-KHI' -2d -auto -opt -portable \
            '--with-unit=physics/Hydro/HydroMain/split/ES' \
            -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -maxblocks=100 \
            -site="$MKFDIR" -objdir="$OBJDIR"
    ;;

    "KHI ES jmark")
        BLOCKSIZE=96
        #MKFDIR=$HOME/projects/turbubox/setups/cheops/flash/makefiles/intel
        MKFDIR=$HOME/projects/turbubox/setups/cheops/flash/makefiles/gcc
        OBJDIR=/mnt/data/flash/KHI/ES/build

        mkdir -p "$OBJDIR"
            
        ./setup 'Markert-KHI' -2d +ug -auto -opt -portable \
            '--with-unit=physics/Hydro/HydroMain/split/ES' \
            -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -maxblocks=-1 \
            -site="$MKFDIR" -objdir="$OBJDIR"
    ;;

    "KHI Bouchut5")
        BLOCKSIZE=96
        MKFDIR=$HOME/turbubox/setups/cheops/flash/makefiles/intel
        OBJDIR=/scratch/jmarker2/flash/KHI/build

        mkdir -p "$OBJDIR"

        ./setup 'Markert-KHI' -2d +ug -auto -opt -portable \
            '--with-unit=physics/Hydro/HydroMain/split/Bouchut5' \
            -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -maxblocks=-1 \
            -site="$MKFDIR" -objdir="$OBJDIR"
    ;;

    "Sedov Blast PPM")
        BLOCKSIZE=96
        MKFDIR="$HOME/turbubox/setups/cheops/flash/makefiles/intel"
        RUNDIR="/scratch/jmarker2/flash/sedov/PPM"
        OBJDIR="${RUNDIR}/build"

        mkdir -p "$RUNDIR/chkpt" "$OBJDIR"

        ./setup 'Sedov' -2d +ug -auto -opt -portable \
            '--with-unit=physics/Hydro/HydroMain/split/PPM' \
            -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -maxblocks=-1 \
            -site="$MKFDIR" -objdir="$OBJDIR"
    ;;

    "Sedov Blast Bouchut5")
        BLOCKSIZE=96
        MKFDIR=$HOME/turbubox/setups/cheops/flash/makefiles/intel
        OBJDIR=/scratch/jmarker2/flash/sedov/build

        mkdir -p "$OBJDIR"

        ./setup 'Sedov' -2d +ug -auto -opt -portable \
            '--with-unit=physics/Hydro/HydroMain/split/Bouchut5' \
            -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -maxblocks=-1 \
            -site="$MKFDIR" -objdir="$OBJDIR"
    ;;

    "Sod Shock b5")
        #BLOCKSIZE=256
        BLOCKSIZE=128
        #BLOCKSIZE=64
        MKFDIR=$HOME/projects/stirturb/turbubox/setup/cheops/makefiles/flash/gcc
        OBJDIR=$HOME/projects/stirturb/simulations/flash/sod-shock/b5/src

        mkdir -p "$OBJDIR"

        ./setup 'Sod' -3d +ug -auto -opt -portable \
            '--with-unit=physics/Hydro/HydroMain/split/Bouchut5' \
            -nxb=$BLOCKSIZE -nyb=4 -nzb=4 -maxblocks=-1 \
            -site="$MKFDIR" -objdir="$OBJDIR"
    ;;

    "Sod Shock b5")
        #BLOCKSIZE=256
        BLOCKSIZE=128
        #BLOCKSIZE=64
        MKFDIR=$HOME/projects/stirturb/turbubox/setup/cheops/makefiles/flash/gcc
        OBJDIR=$HOME/projects/stirturb/simulations/flash/sod-shock/ppm/src

        mkdir -p "$OBJDIR"

        ./setup 'Sod' -3d +ug -auto -opt -portable \
            '--with-unit=physics/Hydro/HydroMain/split/PPM' \
            -nxb=$BLOCKSIZE -nyb=4 -nzb=4 -maxblocks=-1 \
            -site="$MKFDIR" -objdir="$OBJDIR"
    ;;

    "decayturb b3")
        BLOCKSIZE=64
        MKFDIR=$HOME/turbubox/setup/supermuc/makefiles/flash/intel/
        OBJDIR=$HOME/builds/flash/b3/decayturb/

        mkdir -p "$OBJDIR"

        $HOME/turbubox/setup/run-under-intel-env.sh \
            ./setup 'Markert-DecayingTurbulence' -3d +ug -auto -opt -portable \
                '--with-unit=physics/Hydro/HydroMain/split/Bouchut3' \
                -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -nzb=$BLOCKSIZE -maxblocks=-1 \
                -site="$MKFDIR" -objdir="$OBJDIR"
    ;;

    "decayturb b5")
        BLOCKSIZE=64
        MKFDIR=$HOME/turbubox/setup/supermuc/makefiles/flash/intel/
        OBJDIR=$HOME/builds/flash/b5/decayturb/

        mkdir -p "$OBJDIR"

        $HOME/turbubox/setup/run-under-intel-env.sh \
            ./setup 'Markert-DecayingTurbulence' -3d +ug -auto -opt -portable \
                '--with-unit=physics/Hydro/HydroMain/split/Bouchut5' \
                -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -nzb=$BLOCKSIZE -maxblocks=-1 \
                -site="$MKFDIR" -objdir="$OBJDIR"
    ;;

    "stirturb b3")
        MKFDIR=$HOME/turbubox/setup/supermuc/makefiles/flash/intel/
        OBJDIR=$HOME/builds/flash/b3/stirturb/

        mkdir -p "$OBJDIR"

        BLOCKSIZE=64
        $HOME/turbubox/setup/run-under-intel-env.sh \
            ./setup 'Girichidis-StirTurb' -3d +ug -auto -opt -portable \
                '--with-unit=physics/Hydro/HydroMain/split/Bouchut3' \
                -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -nzb=$BLOCKSIZE -maxblocks=-1 \
    ;;

    "stirturb b5")
        BLOCKSIZE=64
        MKFDIR=$HOME/turbubox/setup/supermuc/makefiles/flash/intel/
        OBJDIR=$HOME/builds/flash/b5/stirturb/

        mkdir -p "$OBJDIR"

        $HOME/turbubox/setup/run-under-intel-env.sh \
            ./setup 'Girichidis-StirTurb' -3d +ug -auto -opt -portable \
                '--with-unit=physics/Hydro/HydroMain/split/Bouchut5' \
                -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -nzb=$BLOCKSIZE -maxblocks=-1 \
                -site="$MKFDIR" -objdir="$OBJDIR"
    ;;

    "stirturb ppm")
        BLOCKSIZE=64
        MKFDIR=$HOME/turbubox/setup/supermuc/makefiles/flash/intel/
        OBJDIR=$HOME/builds/flash/ppm/stirturb/

        mkdir -p "$OBJDIR"

        $HOME/turbubox/setup/run-under-intel-env.sh \
            ./setup 'Girichidis-StirTurb' -3d +ug -auto -opt -portable \
                '--with-unit=physics/Hydro/HydroMain/split/PPM' \
                -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -nzb=$BLOCKSIZE -maxblocks=-1 \
                -site="$MKFDIR" -objdir="$OBJDIR"
    ;;

    *)
        echo "Unknown argument: '$arg'" 2>&1
        exit 1

esac; done
