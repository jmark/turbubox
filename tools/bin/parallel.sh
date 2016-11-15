#!/bin/sh

PARALLEL='parallel'

if expr "$(hostname)" : '^cheops' > /dev/null
then
    PARALLEL=$HOME/libs/parallel/src/parallel
fi

$PARALLEL --rpl '{#} 1 $_=sprintf "%04d", $job->seq()-1' "$@"
