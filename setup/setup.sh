#!/usr/bin/env bash

site='cheops.rrz.uni-koeln.de'
#site='cuda.ph1.uni-koeln.de'

NX="${1:?No NX given!}"
solver="${2:?No solver given!}"
objdir="${3:?No source directory given!}"

declare -A solverUnit
solverUnit['8w']='+8wave'
solverUnit['b3']='--with-unit=physics/Hydro/HydroMain/split/Bouchut3'
solverUnit['b5']='--with-unit=physics/Hydro/HydroMain/split/Bouchut5'
solverUnit['es']='--with-unit=physics/Hydro/HydroMain/split/ES'

./setup \
    Girichidis-StirTurb \
    -3d -nxb=$NX -nyb=$NX -nzb=$NX \
    -maxblocks=200 -auto -portable -opt -site="$site" \
    ${solverUnit[$solver]} \
    -objdir="$objdir"
