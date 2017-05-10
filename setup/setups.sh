#!/bin/sh

set -e

function y { true; }
function x { false; }

# mask python3 default installation
if python --version | grep -iq 'python 3'
then
    TMPDIR=$(mktemp -d)
    trap "rm -fr $TMPDIR" EXIT
    ln -sf $(which python2) $TMPDIR/python 
    export PATH="/$TMPDIR:$PATH"
fi

# SOLVER='--with-unit=physics/Hydro/HydroMain/split/PPM'
# SOLVER='--with-unit=physics/Hydro/HydroMain/split/MHD_8Wave'
# SOLVER='--with-unit=physics/Hydro/HydroMain/split/Bouchut3'
# SOLVER='--with-unit=physics/Hydro/HydroMain/split/Bouchut5'
# SOLVER='--with-unit=physics/Hydro/HydroMain/split/ES'

x && { # decayturb b3

    BLOCKSIZE=64
    MKFDIR=$HOME/turbubox/setup/supermuc/makefiles/flash/intel/
    OBJDIR=$HOME/builds/flash/b3/decayturb/

    mkdir -p "$OBJDIR"

    $HOME/turbubox/setup/run-under-intel-env.sh \
        ./setup 'Markert-DecayingTurbulence' -3d +ug -auto -opt -portable \
            '--with-unit=physics/Hydro/HydroMain/split/Bouchut3' \
            -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -nzb=$BLOCKSIZE -maxblocks=-1 \
            -site="$MKFDIR" -objdir="$OBJDIR"
}

x && { # decayturb b5

    BLOCKSIZE=64
    MKFDIR=$HOME/turbubox/setup/supermuc/makefiles/flash/intel/
    OBJDIR=$HOME/builds/flash/b5/decayturb/

    mkdir -p "$OBJDIR"

    $HOME/turbubox/setup/run-under-intel-env.sh \
        ./setup 'Markert-DecayingTurbulence' -3d +ug -auto -opt -portable \
            '--with-unit=physics/Hydro/HydroMain/split/Bouchut5' \
            -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -nzb=$BLOCKSIZE -maxblocks=-1 \
            -site="$MKFDIR" -objdir="$OBJDIR"
}

y && { # stirturb b3

    MKFDIR=$HOME/turbubox/setup/supermuc/makefiles/flash/intel/
    OBJDIR=$HOME/builds/flash/b3/stirturb/

    mkdir -p "$OBJDIR"

    BLOCKSIZE=64
    $HOME/turbubox/setup/run-under-intel-env.sh \
        ./setup 'Girichidis-StirTurb' -3d +ug -auto -opt -portable \
            '--with-unit=physics/Hydro/HydroMain/split/Bouchut3' \
            -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -nzb=$BLOCKSIZE -maxblocks=-1 \
            -site="$MKFDIR" -objdir="$OBJDIR"
}

x && { # stirturb b5

    BLOCKSIZE=64
    MKFDIR=$HOME/turbubox/setup/supermuc/makefiles/flash/intel/
    OBJDIR=$HOME/builds/flash/b5/stirturb/

    mkdir -p "$OBJDIR"

    $HOME/turbubox/setup/run-under-intel-env.sh \
        ./setup 'Girichidis-StirTurb' -3d +ug -auto -opt -portable \
            '--with-unit=physics/Hydro/HydroMain/split/Bouchut5' \
            -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -nzb=$BLOCKSIZE -maxblocks=-1 \
            -site="$MKFDIR" -objdir="$OBJDIR"
}

x && { # stirturb ppm

    BLOCKSIZE=64
    MKFDIR=$HOME/turbubox/setup/supermuc/makefiles/flash/intel/
    OBJDIR=$HOME/builds/flash/ppm/stirturb/

    mkdir -p "$OBJDIR"

    $HOME/turbubox/setup/run-under-intel-env.sh \
        ./setup 'Girichidis-StirTurb' -3d +ug -auto -opt -portable \
            '--with-unit=physics/Hydro/HydroMain/split/PPM' \
            -nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -nzb=$BLOCKSIZE -maxblocks=-1 \
            -site="$MKFDIR" -objdir="$OBJDIR"
}
