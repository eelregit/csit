#!/usr/bin/env python3
# coding: utf-8

import numpy as np
from scipy import integrate
from scipy import interpolate
from scipy.integrate import odeint
from nbodykit.lab import *
from nbodykit import setup_logging

cosmo = cosmology.Planck15
#cosmo_new = cosmo.match(sigma8=0.846687)

seedlist = []
for i in range(1991, 1991+1):
    seedlist.append(str(i))

force = 'treepm2048'
box = 1000

comm = CurrentMPIComm.get()
print(f'rank: {comm.rank} / {comm.size}', flush=True)
setup_logging("debug")

def compute_fof(kind, seed):

    Aniss_data = np.loadtxt('/home/yinli/csit/analysis/Aniss/Aniss_planck2015_%s_z200' % (kind))

    delta_ax = interpolate.interp1d(Aniss_data[:,0], Aniss_data[:,1], kind="cubic")
    delta_ay = interpolate.interp1d(Aniss_data[:,0], Aniss_data[:,2], kind="cubic")
    delta_az = interpolate.interp1d(Aniss_data[:,0], Aniss_data[:,3], kind="cubic")

    def a_ratio(z):
        scale_a = 1./(1+z)
        # replace following with function of redshift
        delta_a = np.zeros(3)
        delta_a[0] = delta_ax(scale_a)
        delta_a[1] = delta_ay(scale_a)
        delta_a[2] = delta_az(scale_a)
        return 1 + delta_a

    if comm.rank == 0:
        print(comm.size)
        print(kind)
        print('reading snapshot files...')
#    part_cat = Gadget1Catalog('/mnt/sdceph/users/yinli/csit/planck2015/%d/np1024/%s/%s/nbody/treepm2048/snapdir_005/snap_005.*' % (box, seed, kind))
    part_cat = Gadget1Catalog('/mnt/sdceph/users/yinli/csit/planck2015/%d/%s/%s/SIMP_ASMTH45/nbody/snapdir_003/snap_003.*' % (box, seed, kind))

    if comm.rank == 0:
        print(part_cat)
        print(part_cat.columns)
    # rescale coordinates
    z = part_cat.attrs['Redshift']
    part_cat.attrs['BoxSize'] *= a_ratio(z)  # from kpc/h to Mpc/h
    part_cat['Position'] *= a_ratio(z)
    part_cat['Velocity'] = part_cat['GadgetVelocity'] * a_ratio(z) ** 0.5  ##if Velocity is peculiar velocity,then v_p = a*dx/dt and GadgetVelocity u = v_p/sqrt(a).

    part_cat.attrs['Nmesh'] = (1,) * part_cat['Position'].shape[1]  # HACK to make FOF compute the mean_separation

    part_mass = part_cat['Mass'][0] * 1e10  # Msun/h

#    f = FOF(part_cat, 0.05, 20, absolute=True)
    f = FOF(part_cat, 0.2, 20, absolute=True)

    halo_cat = f.to_halos(part_mass, cosmo, 0)  # NOTE check if this is the right cosmology

    if comm.rank == 0:
        print(halo_cat)
        print(halo_cat.columns)

#    halo_cat.save('/mnt/sdceph/users/yinli/csit/planck2015/%d/np1024/%s/%s/halos' % (box, seed, kind))
    halo_cat.save('/mnt/sdceph/users/yinli/csit/planck2015/%d/%s/%s/SIMP_ASMTH45/halos' % (box, seed, kind))


if __name__ == "__main__":
    for seed_num in seedlist:
#        for filename in ['iso', 'mmp001', 'ppm001', 'mmp002', 'ppm002', 'mmp005', 'ppm005', 'mmp008', 'ppm008' , 'mmp01', 'ppm01', 'mmp02', 'ppm02', 'mmp03', 'ppm03' , 'mmp04', 'ppm04', 'mmm0004', 'ppp0004', 'mmm0008', 'ppp0008', 'mmm001', 'ppp001','mmm002', 'ppp002','mmm003', 'ppp003','mmm007', 'ppp007','mmm01', 'ppp01']:
#        for filename in ['iso', 'ppp0004', 'mmm0004', 'ppp007', 'mmm007', 'mmp01', 'ppm01', 'mmp001', 'ppm001']:
#       for filename in ['mmm003', 'ppp007', 'mmm007', 'ppp0004', 'mmm0004', 'mmp01', 'ppm01', 'm01m01p01', 'p01p01m01', 'm002m002p013', 'm008m008p007', ]:
        for filename in ['iso', 'ppp003', 'mmm003', 'ppp007', 'mmm007', 'ppp0004', 'mmm0004']:
            compute_fof(filename, seed_num)
