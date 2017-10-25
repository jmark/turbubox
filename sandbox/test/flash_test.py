#!/usr/bin/env python3

import flash
import sys

fname = sys.argv[1]
fl = flash.File(fname)
print(fl.gridsize)
