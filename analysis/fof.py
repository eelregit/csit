#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from nbodykit.lab import *

from nbodykit import setup_logging, style
from nbodykit import use_mpi
from nbodykit import CurrentMPIComm
from nbodykit.io.gadget import Gadget1File
from nbodykit.lab import ArrayMesh
from nbodykit.source.catalog import Gadget1Catalog

import math
from scipy import integrate
from scipy import interpolate
from scipy.integrate import odeint

cosmo = cosmology.Planck15

seed = 2006
force = 'treepm2048'
box = 1000

comm = CurrentMPIComm.get()
use_mpi(comm=comm)
setup_logging("debug")

Omega_m = cosmo.Omega_m(0)
Omega_L = cosmo.Omega_lambda(0)
tau_array = [0.005, 0.005, -0.01]
a_begin = 1/(1+200)
a_end = 1
points = 5000

def Hubble(a):
    matter = Omega_m/a**3
    Lambda = Omega_L
    return np.sqrt(matter+Lambda)

def dHda(a):
    up = -3*Omega_m
    down = 2*a**4*Hubble(a)
    return up/down

def dHdloga(a):
    return a*dHda(a)

def growth_int(a):
    return 1./(a*Hubble(a))**3

def growth_D(a):
    factor = 5*Omega_m*Hubble(a)/2
    return factor*integrate.quad(growth_int, 0, a)[0]

def dDda(a):
    factor = 5*Omega_m*dHda(a)/2
    second = 5*Omega_m/(2*a**3*Hubble(a)**2)
    return factor*integrate.quad(growth_int, 0, a)[0] + second

def dDdloga(a):
    return a*dDda(a)

vgrowth_D = np.vectorize(growth_D)
vdDda = np.vectorize(dDda)
vdDdloga = np.vectorize(dDdloga)

a_array = np.logspace(np.log(a_begin), np.log(1), points, base=np.exp(1))
lna_array = np.log(a_array)
a_da = np.zeros((points,2,3))

def Delta_f(x, lna, tau):
    a = np.exp(lna)
    first = 2+dHdloga(a)/Hubble(a)
    timetau = tau*vgrowth_D(a)/vgrowth_D(1)
    const = 3*Omega_m*timetau/(2*a**3*Hubble(a)**2)
    ret = [
        x[1],
        -first*x[1] -const*(1+x[0])
    ]
    return ret

inifac = growth_D(a_array[0])/growth_D(1)
dinifac = inifac*vdDdloga(a_array[0])/growth_D(a_array[0])

inifac_array = -np.array([inifac, dinifac])

Delta_ini = np.array([tau_array[0]*inifac_array,
                      tau_array[1]*inifac_array,
                      tau_array[2]*inifac_array])

for i in range(3):
    a_da[:,:,i] = odeint(Delta_f, Delta_ini[i], lna_array, args=(tau_array[i],))

delta_a_func_x = interpolate.interp1d(a_array, a_da[:,0,0], kind="cubic")
delta_a_func_y = interpolate.interp1d(a_array, a_da[:,0,1], kind="cubic")
delta_a_func_z = interpolate.interp1d(a_array, a_da[:,0,2], kind="cubic")

def a_ratio_mmp(z):
    scale_a = 1./(1+z)
    # replace following with function of redshift
    delta_a = np.zeros(3)
    delta_a[0] = delta_a_func_x(scale_a)
    delta_a[1] = delta_a_func_y(scale_a)
    delta_a[2] = delta_a_func_z(scale_a)
    return 1 + delta_a

def a_ratio_ppm(z):
    scale_a = 1./(1+z)
    # replace following with function of redshift
    delta_a = np.zeros(3)
    delta_a[0] = delta_a_func_x(scale_a)
    delta_a[1] = delta_a_func_y(scale_a)
    delta_a[2] = delta_a_func_z(scale_a)
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
    if kind=='mmp':
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

for filename in ['iso', 'mmp', 'ppm']:
    compute_fof(filename)
