#!/usr/bin/env python3
# coding: utf-8

import numpy as np
from nbodykit.lab import *
from nbodykit import setup_logging
from nbodykit.source.catalog import Gadget1Catalog

box = 1000
seed = 1991

comm = CurrentMPIComm.get()
print(f'rank: {comm.rank} / {comm.size}', flush=True)
setup_logging()

def compute_power(kind):
    num_mesh = 2048
    if comm.rank == 0:
        print(comm.size)
        print(kind)
        print('reading snapshot files...')
    part_cat = Gadget1Catalog('/mnt/sdceph/users/yinli/csit/ic/planck2015/%d/%d/%s/ic.*' % (box, seed, kind))
    if comm.rank == 0:
        print(part_cat)
        print(part_cat.columns)
    mesh = part_cat.to_mesh(resampler='pcs', Nmesh=num_mesh, BoxSize=box, compensated=True, position='Position', interlaced=True)
    if comm.rank == 0:
        print('starting FFT...')
    FFT = FFTPower(mesh, mode='2d', dk=0.04, Nmu=200, los=[0,0,1], poles=[0,2,4])
    if comm.rank == 0:
        print('done FFT')
    FFT.save("/home/yinli/csit/analysis/power/planck2015/%d/%d/ic_%s_fftpower.json" % (box, seed, kind))
    poles = FFT.poles
    mono = poles['power_0'].real
    quad = poles['power_2'].real
    hexa = poles['power_4'].real
    p0 = np.array([poles['k'][1:], mono[1:]])
    p2 = np.array([poles['k'][1:], quad[1:]])
    p4 = np.array([poles['k'][1:], hexa[1:]])
    np.savetxt('/home/yinli/csit/analysis/power/planck2015/%d/%d/ic_%s_p0.dat' % (box, seed, kind), p0)
    np.savetxt('/home/yinli/csit/analysis/power/planck2015/%d/%d/ic_%s_p2.dat' % (box, seed, kind), p2)
    np.savetxt('/home/yinli/csit/analysis/power/planck2015/%d/%d/ic_%s_p4.dat' % (box, seed, kind), p4)

if __name__ == "__main__":
    for filename in ['iso', 'mmp', 'ppm']:
        compute_power(filename)
