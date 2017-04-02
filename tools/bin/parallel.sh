#!/bin/sh

PARALLEL='parallel'

if hostname | grep -qE '^cheops'
then
    PARALLEL=$HOME/libs/parallel/src/parallel
fi

$PARALLEL --rpl '{#} 1 $_=sprintf "%04d", $job->seq()-1' "$@"
