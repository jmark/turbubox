#!/usr/bin/env bash

set -euo pipefail

function workdir
{
    local grid=$1; shift
    local unit=$1; shift
    local solv=$1; shift

    echo -n "/mnt/seagate/FLASH/stirturb/mach-$MACH/$unit/$grid/$solv"
}
export -f workdir

function setup
{
    local grid="$1"; shift
    local unit="$1"; shift
    local solv="$1"; shift
    local wkdr="$(workdir $grid $unit $solv)"

    mkdir -vp $wkdr/{src,amr,run,log,csv}

    local refine=2
    local NX=$(expr $grid / $refine)

    declare -A solverUnit
    solverUnit['8w']='+8wave'
    solverUnit['b3']='--with-unit=physics/Hydro/HydroMain/split/Bouchut3'
    solverUnit['b5']='--with-unit=physics/Hydro/HydroMain/split/Bouchut5'
    solverUnit['es']='--with-unit=physics/Hydro/HydroMain/split/ES'

    pushd '/srv/projects/astro/frameworks/silcc/code'
        ./setup \
            Girichidis-StirTurb \
            -3d -nxb=$NX -nyb=$NX -nzb=$NX \
            -maxblocks=200 -auto -opt -site=archlinux-intel \
            -portable \
            ${solverUnit[$solv]} \
            -objdir="$wkdr/src"

        log_exit_status ${grid} ${unit} ${solv} $?
    popd
}
export -f setup

function compile
{
    local grid="$1"; shift
    local unit="$1"; shift
    local solv="$1"; shift
    local wkdr="$(workdir $grid $unit $solv)"
    
    pushd "$wkdr/src" 
        make -j12

        log_exit_status ${grid} ${unit} ${solv} $?
    popd
}
export -f compile

function runtime
{
    local grid="$1"; shift
    local unit="$1"; shift
    local solv="$1"; shift
    local wkdr="$(workdir $grid $unit $solv)"

    pushd "$wkdr/run" 
        ln -vfrs ../src/flash4 .
    popd

    ./mk-config.$unit.tcl flash.par.ttcl $MACH -json \
         > "$wkdr/run/flash.par" \
        2> "$wkdr/run/META.json"

    log_exit_status ${grid} ${unit} ${solv} $?
}
export -f runtime 

function run
{
    local grid="$1"; shift
    local unit="$1"; shift
    local solv="$1"; shift
    local wkdr="$(workdir $grid $unit $solv)"

    local fmt='E e S U P M t K D p X F R W c w k Z x'

    echo "Running $grid/$unit/$solv..."
    echo "tail -f $wkdr/run/stdout+stderr.log" > /tmp/tail-f.sh
    echo "Follow process via: bash /tmp/tail-f.sh"

    pushd "$wkdr/run" > /dev/null
        set +e
        /usr/bin/time -f %${fmt// /\\t%} -o TIMING.DAT \
            /opt/mpich/bin/mpirun -n $(nproc) ./flash4 > stdout+stderr.log 

        log_exit_status ${grid} ${unit} ${solv} $?
        set -e
    popd > /dev/null

    echo
}
export -f run

function analyze
{
    local grid="$1"; shift
    local unit="$1"; shift
    local solv="$1"; shift
    local wkdr="$(workdir $grid $unit $solv)"

    pushd "$wkdr/run"
        CSV=$(cat "TIMING.DAT")
    popd

    echo "grid\tunit\t\solv\t$CSV" >> /tmp/timing.dat

    # ...
}
export -f analyze

function log_exit_status
{
    local grid="$1"; shift
    local unit="$1"; shift
    local solv="$1"; shift
    local stat="$1"; shift

    echo -e "$grid\t$unit\t$solv\t$stat" | tee -a $LOGFILE
}
export -f log_exit_status

export SCRIPTNAME="$0"
export COMMAND="$1"
export LOGFILE="/tmp/$SCRIPTNAME.log"

echo >> $LOGFILE 
echo "${COMMAND}: $(date --rfc-3339=seconds)" >> $LOGFILE

export MACH=20.0

grids=(16 24 32 48 64)
units=(unit cgs)
solvs=(8w b3 b5 es)

grids=(16)
units=(unit)
solvs=(es)

parallel -j1 --ungroup \
    $COMMAND ::: ${grids[@]} ::: ${units[@]} ::: ${solvs[@]} 
