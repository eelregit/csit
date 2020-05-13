#!/usr/bin/env python3
# coding: utf-8

import numpy as np
from nbodykit.lab import *
from nbodykit import setup_logging
from nbodykit.source.catalog import Gadget1Catalog

force = 'treepm2048'
box = 1000
seedlist = []
for i in range(1991,1991+1):
    seedlist.append(str(i))

comm = CurrentMPIComm.get()
print(f'rank: {comm.rank} / {comm.size}', flush=True)
setup_logging()

def dh(kind):
    if kind=='iso':
        return 1.
    elif kind=='ppp0004':
        return 0.674727/0.6774
    elif kind=='mmm0004':
        return 0.680063/0.6774
    elif kind=='ppp003':
        return 0.657085/0.6774
    elif kind=='mmm003':
        return 0.697123/0.6774
    elif kind=='ppp007':
        return 0.628978/0.6774
    elif kind=='mmm007':
        return 0.722584/0.6774

def compute_power(kind, seed):
    num_mesh = 8196
    if comm.rank == 0:
        print(comm.size)
        print(kind)
        print('reading snapshot files...')
#    part_cat = Gadget1Catalog('/mnt/sdceph/users/yinli/csit/planck2015/%d/np1024/%s/%s/nbody/treepm2048/snapdir_005/snap_005.*' % (box, seed, kind))
    part_cat = Gadget1Catalog('/mnt/sdceph/users/yinli/csit/planck2015/%d/%s/%s/SU/nbody/snapdir_005/snap_005.*' % (box, seed, kind))

    if comm.rank == 0:
        print(part_cat)
        print(part_cat.columns)
    mesh = part_cat.to_mesh(resampler='pcs', Nmesh=num_mesh, BoxSize=box, compensated=True, position='Position', interlaced=True)
    if comm.rank == 0:
        print('starting FFT...')
    FFT = FFTPower(mesh, mode='2d', dk=0.015/dh(kind), Nmu=200, los=[0,0,1], poles=[0,2,4])
    if comm.rank == 0:
        print('done FFT')
#    FFT.save('/mnt/sdceph/users/yinli/csit/planck2015/%d/np1024/%s/%s/power/z0_pm4096_fftpower.json' % (box, seed, kind))
    FFT.save('/mnt/sdceph/users/yinli/csit/planck2015/%d/%s/%s/SU/power/z0_fftpower.json' % (box, seed, kind))
    poles = FFT.poles
    mono = poles['power_0'].real
    quad = poles['power_2'].real
    hexa = poles['power_4'].real
    p0 = np.array([poles['k'][1:], mono[1:]])
    p2 = np.array([poles['k'][1:], quad[1:]])
    p4 = np.array([poles['k'][1:], hexa[1:]])
#    np.savetxt('/mnt/sdceph/users/yinli/csit/planck2015/%d/np1024/%s/%s/power/z0_treepm2048_p0.dat' % (box, seed, kind), p0)
    np.savetxt('/mnt/sdceph/users/yinli/csit/planck2015/%d/%s/%s/SU/power/z0_p0.dat' % (box, seed, kind), p0)
    np.savetxt('/mnt/sdceph/users/yinli/csit/planck2015/%d/%s/%s/SU/power/z0_p2.dat' % (box, seed, kind), p2)
    np.savetxt('/mnt/sdceph/users/yinli/csit/planck2015/%d/%s/%s/SU/power/z0_p4.dat' % (box, seed, kind), p4)

if __name__ == "__main__":
    for seed_num in seedlist:
#        for filename in ['mmp002', 'ppm002', 'mmp005', 'ppm005', 'mmp008', 'ppm008', 'mmp02', 'ppm02', 'mmp03', 'ppm03' , 'mmp04', 'ppm04','mmm0008', 'ppp0008', 'mmm001', 'ppp001','mmm002', 'ppp002','mmm003', 'ppp003','mmm01', 'ppp01']:
#       for filename in ['iso', 'ppp003', 'mmm003', 'ppp007', 'mmm007', 'ppp0004', 'mmm0004']:
       for filename in ['iso']:
            compute_power(filename, seed_num)
