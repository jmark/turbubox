#!/usr/bin/env python3

import numpy as np
import sys

import ulz
import flexi
import hopr

hoprfp = sys.argv[1]

for count,flxfp in enumerate(sys.stdin):
    flx = flexi.File(flxfp.rstrip(), hopr.CartesianMeshFile(hoprfp))

    time = flx.time

    # TODO
    rho0  = 1.0 
    kappa = 5/3
    mu0   = 1.0

    dens, velx, vely, velz, pres = flx.get_prims()
    
    c_s = np.sqrt(np.mean(pres)/np.mean(dens))

    GS = np.array(dens.shape)
    DS = flx.hopr.domainsize
    CS = DS/GS

    Vgrid   = np.prod(GS)
    Vcell   = np.prod(CS) 
    Vdomain = np.prod(DS) 

    vortx,vorty,vortz = ulz.curl(velx,vely,velz, CS[0], CS[1], CS[2])

    masstotal = np.sum(dens) * Vdomain
    ekintotal = Vcell/2.0 * np.sum(dens * (velx**2+vely**2+velz**2))
    einttotal = np.mean(pres)/(kappa-1) 
    etotal    = ekintotal + einttotal

    machglob  = np.sqrt(np.sum(dens * (velx**2+vely**2+velz**2))/np.sum(dens)) / c_s

    evortotal = CS[0]**5/12.0 * np.sum(dens * (vortx**2+vorty**2+vortz**2))

    print("\t".join(map(str,[
        count, time, mach etotal, einttotal, ekintotal, evortotal
    ])), flush=True)
