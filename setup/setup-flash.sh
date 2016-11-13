#!/bin/sh

set -eu

if expr "$1" : 'help' > /dev/null
then
	echo "usage: setup.sh <block size> <max blocks> <solver> <objdir> <site dir> <additional options>"
	exit 1
fi

BLOCKSIZE="${1:?No block size given!}"          && shift
MAXBLOCKS="${1:?No maxblocks given}"            && shift
   SOLVER="${1:?No solver given!}"              && shift
   OBJDIR="${1:?No source directory given!}"    && shift
  SITEDIR="${1:?No site directory given!}"      && shift

declare -A SOLVERUNIT
SOLVERUNIT['pm']='--with-unit=physics/Hydro/HydroMain/split/PPM'
SOLVERUNIT['8w']='--with-unit=physics/Hydro/HydroMain/split/MHD_8Wave'
SOLVERUNIT['b3']='--with-unit=physics/Hydro/HydroMain/split/Bouchut3'
SOLVERUNIT['b5']='--with-unit=physics/Hydro/HydroMain/split/Bouchut5'
SOLVERUNIT['es']='--with-unit=physics/Hydro/HydroMain/split/ES'

# mask python3 default installation
if expr "$(hostname)" : '^jmark' > /dev/null
then
    TMPDIR=$(mktemp -d)
    trap "rm -r $TMPDIR" EXIT
    ln -sf $(which python2) $TMPDIR/python 
    export PATH="/$TMPDIR:$PATH"
fi

./setup Girichidis-StirTurb \
    -3d -auto -portable -opt \
	-nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -nzb=$BLOCKSIZE -maxblocks=$MAXBLOCKS \
	-site="$SITEDIR" -objdir="$OBJDIR" ${SOLVERUNIT[$SOLVER]} \
    "$@"
