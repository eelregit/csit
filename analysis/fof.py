#!/usr/bin/env python3
# coding: utf-8

import numpy as np
from scipy import integrate
from scipy import interpolate
from scipy.integrate import odeint
from nbodykit.lab import *
from nbodykit import setup_logging

cosmo = cosmology.Planck15

seed = 2006
force = 'treepm2048'
box = 1000

comm = CurrentMPIComm.get()
print(f'rank: {comm.rank} / {comm.size}', flush=True)
setup_logging("debug")

data_mmp = np.loadtxt('/home/yinli/csit/lgadget2/Aniss/Aniss_planck2015_m0005m0005p001_z200')
#data_ppm = np.loadtxt('/home/yinli/csit/lgadget2/Aniss/Aniss_planck2015_p0005p0005m001_z200')

delta_ax = interpolate.interp1d(data_mmp[0], data_mmp[1], kind="cubic")
delta_ay = interpolate.interp1d(data_mmp[0], data_mmp[2], kind="cubic")
delta_az = interpolate.interp1d(data_mmp[0], data_mmp[3], kind="cubic")

def a_ratio_mmp(z):
    scale_a = 1./(1+z)
    # replace following with function of redshift
    delta_a = np.zeros(3)
    delta_a[0] = delta_ax(scale_a)
    delta_a[1] = delta_ay(scale_a)
    delta_a[2] = delta_az(scale_a)
    return 1 + delta_a

def a_ratio_ppm(z):
    scale_a = 1./(1+z)
    # replace following with function of redshift
    delta_a = np.zeros(3)
    delta_a[0] = delta_ax(scale_a)
    delta_a[1] = delta_ay(scale_a)
    delta_a[2] = delta_az(scale_a)
    return 1 - delta_a

def compute_fof(kind):
    if comm.rank == 0:
        print(comm.size)
        print(kind)
        print('reading snapshot files...')

    part_cat = Gadget1Catalog('/mnt/sdceph/users/yinli/csit/nbody/planck2015/%d/%s/%d/%s/snapdir_004/snap_004.*' % (box, force, seed, kind))

    if comm.rank == 0:
        print(part_cat)
        print(part_cat.columns)
    # rescale coordinates
    z = part_cat.attrs['Redshift']
    if kind=='iso':
        part_cat['Velocity'] = part_cat['GadgetVelocity'] * np.ones(3)
    elif kind=='mmp':
        part_cat.attrs['BoxSize'] /= a_ratio_mmp(z)  # from kpc/h to Mpc/h
        part_cat['Position'] /= a_ratio_mmp(z)
        part_cat['Velocity'] = part_cat['GadgetVelocity'] * a_ratio_mmp(z) ** 0.5  ##Velocity is peculiar velocity, i.e. v_p = a*dx/dt
    elif kind=='ppm':
        part_cat.attrs['BoxSize'] /= a_ratio_ppm(z)
        part_cat['Position'] /= a_ratio_ppm(z)
        part_cat['Velocity'] = part_cat['GadgetVelocity'] * a_ratio_ppm(z) ** 0.5  ##Velocity is peculiar velocity, i.e. v_p = a*dx/dt

    part_cat.attrs['Nmesh'] = (1,) * part_cat['Position'].shape[1]  # HACK to make FOF compute the mean_separation

    part_mass = part_cat['Mass'][0] * 1e10  # Msun/h

    f = FOF(part_cat, 0.2, 20)

    halo_cat = f.to_halos(part_mass, cosmo, 0)  # NOTE check if this is the right cosmology

    if comm.rank == 0:
        print(halo_cat)
        print(halo_cat.columns)

    halo_cat.save('home/yinli/csit/analysis/halos/planck2015/%d/%d/halo_%s_%s' % (box, seed, force, kind))

if __name__ == "__main__":
    for filename in ['iso', 'mmp', 'ppm']:
        compute_fof(filename)
