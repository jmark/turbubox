#!/usr/bin/env bash

NX=16

solv='es'
site='archlinux-intel'

wkdr="$(pwd)/.."

declare -A solverUnit
solverUnit['8w']='+8wave'
solverUnit['b3']='--with-unit=physics/Hydro/HydroMain/split/Bouchut3'
solverUnit['b5']='--with-unit=physics/Hydro/HydroMain/split/Bouchut5'
solverUnit['es']='--with-unit=physics/Hydro/HydroMain/split/ES'

pushd '/srv/projects/astro/frameworks/silcc/code'
    ./setup \
        Girichidis-StirTurb \
        -3d -nxb=$NX -nyb=$NX -nzb=$NX \
        -maxblocks=200 -auto -portable -opt -site=$site \
        ${solverUnit[$solv]} \
        -objdir="$wkdr/src"
popd
