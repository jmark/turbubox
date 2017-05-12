#!/bin/bash

set -u

function y { true; }
function x { false; }

snapshots() {
    local FILES
    FILES=""
    FILES="$FILES $(test -e checkpoints && find checkpoints/ -name 'flash_hdf5_plt_cnt_*' | sort)"
    FILES="$FILES $(find . -name 'sim_State_*.h5' | sort)"
    echo $FILES
}

DIRS="$@"
CONF="$HOME/turbubox/setup/supermuc/loadleveler/hw.micro.conf"

for DIR in $DIRS
do
pushd $DIR

    SETUP=$(basename $DIR)
    
    x && {
        rm  -fv *.log *.sh *.ini *.tmp *.conf *.TXT
        rm -rfv pickle*
        rm -rfv png*
        rm -rfv cache*
        rm -rfv tmp log checkpoints
    }

    x && {
        ~/scripts/parallel/llrun.pl -f ../sb.general.flexi.conf \
            :: --job_name "flexi-sim: stirturb/$SETUP" \
            :: poe $HOME/builds/flexi/builds/ibm-mpi/bin/flexi flexi.ini \
            | llsubmit -
    }

    x && { ## analysis
        mkdir -p log pickle
        FILES=$(snapshots)

        ~/scripts/parallel/llrun.pl -f "$CONF" \
            :: --job_name "anal-proc: stirturb/$SETUP" --output "log/ana.log" --error "log/ana.log" \
               --wall_clock_limit 12:00:00 \
            :: $HOME/turbubox/tools/bin/analysis.py --destdir pickle/ --skip --parallel ${NPROCS:-2} $FILES \
            | llsubmit -
    }
    
    x && { ## plotting

        CACHEDIR='cache'
        PNGDIR='png'
        LOGDIR='log'

        mkdir -p $LOGDIR $PNGDIR $CACHEDIR

        FILES=$(snapshots)
        NFILES="$(echo "$FILES" | wc -w)"
        #CRANGE="cdens=(-1.0,1.0), cekin=(0.02,2.0), cmach=(0.2,1.2), cvort=(-12,-7)"
        CRANGE="cdens=(-1.0,1.0), cekin=(0.02,1.0), cmach=(0.2,1.0), cvort=(-10,-7)"

        ~/scripts/parallel/llrun.pl -f "$CONF" \
            :: --job_name "plot-proc: stirturb/$SETUP" --output "$LOGDIR/plot.log" --error "$LOGDIR/plot.log" \
               --wall_clock_limit 12:00:00 \
            :: $HOME/turbubox/plot/dens-ekin-mach-vort.py --destdir "$PNGDIR" --cachedir "$CACHEDIR" \
                --parallel ${NPROCS:-2} --ntasks $NFILES --crange "$CRANGE" --title "$SETUP" \
                $FILES | llsubmit -
    }

    y && { # progress counting
        FILES="$(test -e checkpoints && find checkpoints/ -name 'flash_hdf5_plt_cnt_*' | sort)"
        FILES="$FILES $(find . -name 'sim_State_*.h5' | sort)"
        NFILES="$(echo "$FILES" | wc -w)"
        echo    '#snapshots:' $NFILES
        echo -n '#pickle:    ' && find pickle/ -name '*.pickle' | wc -l
        echo -n '#png:       ' && find png/ -name '*.png' | wc -l
        echo -n '#tmp:       ' && find . -name '*.tmp' | wc -l
    }

    x && {
        grep 'st_mach' flexi.ini
        grep 'bulkmotion%' flexi.ini
    }

popd > /dev/null
echo
done
