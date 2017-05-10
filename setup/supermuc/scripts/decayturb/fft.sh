#!/bin/bash

mkdir -p pickle

$HOME/turbubox/tools/bin/analysis.py --destdir pickle/ --skip --parallel ${NPROCS:-4} "$@"
