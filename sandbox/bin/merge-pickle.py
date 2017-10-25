#!/usr/bin/env python3

# stdlib
import os
import sys
import pickle
import numpy as np

# jmark
import dslopts

with dslopts.Manager(scope=globals(),appendix="pickle files are be defined after '--'.") as mgr:
    mgr.add('sinkfp')

rows = []
for fp in ARGV_TAIL:
    with open(fp, 'rb') as fd:
        row = pickle.load(fd)
        rows.append(row)
    print(fp)

result = rows

with open(sinkfp, 'wb') as fd:
    pickle.dump(result, fd)
