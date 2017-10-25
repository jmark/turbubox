#!/bin/sh

exec ./gcc-env.sh ./build/couchdg \
    --destdir   data \
    --analdir   . \
    --inittime  0.0 \
    --stoptime  0.5 \
    --dtchkpt   0.01 \
    #--initfile  data/chkpt_0126.h5
