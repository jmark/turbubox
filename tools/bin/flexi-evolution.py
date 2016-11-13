#!/usr/bin/env python3

import numpy as np
import sys

import flash as fl
import ulz

import flexi
import hopr

hoprfp = sys.argv[1]

for count,flxfp in enumerate(sys.stdin):
    flx = flexi.File(flxfp.rstrip(), hopr.CartesianMeshFile(hoprfp))

    time = flx.time

    # TODO
    rho0 = 1.0 
    kappa = flx.params['kappa']
    mu0   = 1.0

    cons  = [flx.flexi_to_box(i) for i in range(0,8)]
    prims = ulz.mhd_conservative_to_primitive(cons, kappa, mu0)

    dens  = prims[0] 
    vels  = prims[1:4] 
    pres  = prims[4]
    
    c_s   = np.sqrt(np.mean(pres)/np.mean(dens))

    GS = np.array(dens.shape)
    DS = flx.hopr.domainsize
    CS = DS/GS

    Vgrid   = np.prod(GS)
    Vcell   = np.prod(CS) 
    Vdomain = np.prod(DS) 

    vorts = ulz.curl(*vels,*tuple(CS))

    masstotal = np.sum(dens) * Vdomain
    ekintotal = Vcell/2.0 * np.sum(dens * ulz.norm(*vels))

    mach      = np.sqrt(np.sum(dens * ulz.norm(*vels))/np.sum(dens)) / c_s
    etotal    = np.sum(cons[4])

    evortotal = CS[0]**5/12.0 * np.sum(dens * ulz.norm(*vorts))

    print("\t".join(map(str,[

        count,time, mach, etotal, ekintotal, evortotal, 
        np.mean(prims[4]), np.min(prims[4]), np.max(prims[4])

    ])), flush=True)
