#!/usr/bin/env python

import numpy as np
import interpolate as itpl

input = np.ones((8,8))

output = itpl.foo(input)

print(output)
