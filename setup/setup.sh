#!/usr/bin/env bash

if test "$1" = "help"
then
	echo "usage: setup.sh <block size> <max blocks> <solver> <objdir> <site dir>"
	exit 1
fi

BLOCKSIZE="${1:?No block size given!}"
MAXBLOCKS="${2:?No maxblocks given}"
SOLVER="${3:?No solver given!}"
OBJDIR="${4:?No source directory given!}"
SITEDIR="${5:?No site directory given!}"

declare -A SOLVERUNIT
SOLVERUNIT['8w']='--with-unit=physics/Hydro/HydroMain/split/MHD_8Wave'
SOLVERUNIT['b3']='--with-unit=physics/Hydro/HydroMain/split/Bouchut3'
SOLVERUNIT['b5']='--with-unit=physics/Hydro/HydroMain/split/Bouchut5'
SOLVERUNIT['es']='--with-unit=physics/Hydro/HydroMain/split/ES'

./setup Girichidis-StirTurb -3d -auto -portable -opt \
	-nxb=$BLOCKSIZE -nyb=$BLOCKSIZE -nzb=$BLOCKSIZE -maxblocks=$MAXBLOCKS \
	-site="$SITEDIR" -objdir="$OBJDIR" ${SOLVERUNIT[$SOLVER]}
