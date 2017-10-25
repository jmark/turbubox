#!/bin/bash

function y { true; }
function x { false; }

set -u

DIRS="$@"
CONF="$HOME/turbubox/setup/supermuc/loadleveler/hw.micro.conf"

snapshots() {
    local FILES
    FILES="$(test -e checkpoints && echo checkpoints/flash_forced_hdf5_plt_cnt_0000)"
    FILES="$FILES $(test -e checkpoints && find checkpoints/ -name 'flash_hdf5_plt_cnt_*' | sort)"
    FILES="$FILES $(find . -name 'sim_State_*.h5' | sort)"
    echo $FILES
}

for dir in $DIRS
do
pushd $dir

    SETUP=$(basename $dir)

    y && {
        mkdir -p ~/projects/stirturb/plots/decayturb/snapshots/files/$SETUP
        cp -v png/*.png ~/projects/stirturb/plots/decayturb/snapshots/files/$SETUP
    }

    x && {
        cp -v ../../stirturb/$dir/*.ini .
        cp -v ../rk3-fv/hopr_mesh.h5 .
        cp -v ../rk3-fv/sim_State_0000000.000000000.h5 .
    }

    x && {
        diff ../../stirturb/$dir/flexi.ini flexi.ini
    }

    x && {
        $HOME/scripts/parallel/llrun.pl -f ../sb.general.flexi.conf \
            :: \
            :: poe $HOME/flexi/builds/mpi-ibm/bin/flexi flexi.ini sim_State_0000000.000000000.h5 \
        | llsubmit -
    }

    x && rm -fv pickle/*.pickle
    x && rm -fv pickle/*.tmp
    x && rm -fv png/*.png
    
    x && { ## analysis
        mkdir -p log pickle

        FILES=$(snapshots)

        ~/scripts/parallel/llrun.pl -f "$CONF" \
            :: --job_name "PWS post-processing: $dir" --output "log/ana.log" --error "log/ana.log" \
               --wall_clock_limit 08:00:00 \
            :: $HOME/turbubox/tools/bin/analysis.py --destdir pickle/ --skip --parallel ${NPROCS:-2} $FILES \
            | llsubmit -
    }

    x && { ## plotting
        mkdir -p log cache png

        FILES=$(snapshots)
        NFILES="$(echo "$FILES" | wc -w)"
        CRANGE="cdens=(-1.0,1.0), cekin=(0.02,2.0), cmach=(0.2,1.2), cvort=(-12,-7)"

        ~/scripts/parallel/llrun.pl -f "$CONF" \
            :: --job_name "plotting: $dir" --output "log/plot.log" --error "log/plot.log" \
               --wall_clock_limit 08:00:00 \
            :: $HOME/turbubox/plot/dens-ekin-mach-vort.py --destdir png/ --cachedir cache/ \
                --parallel ${NPROCS:-2} --ntasks $NFILES --crange "$CRANGE" --title "turbulent decay: $dir" \
                $FILES | llsubmit -
    }

    x && { # process counting
        #echo -n '#png:    '
        #find png -name '*.png' | wc -l
        echo -n '#pickle: '
        find pickle -name '*.pickle' | wc -l
    }

    x && { # process counting

        FILES=$(snapshots)
        NFILES="$(echo "$FILES" | wc -w)"

        echo    '#snapshots:' $NFILES
        echo -n '#pickle:    ' && find pickle/ -name '*.pickle' | wc -l
        echo -n '#tmp:       ' && find pickle/ -name '*.tmp' | wc -l
        echo -n '#cache:     ' && find cache/ -name '*.pickle' | wc -l
        echo -n '#png:       ' && find png/ -name '*.png' | wc -l
    }


popd > /dev/null
echo
done
