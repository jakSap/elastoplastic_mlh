#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import h5py

"""
this program creates two 2-D rings (z = 0) around the origin which are then shifted to their final positions
used for colliding rings testcase, see Gray, Monaghan, Swift SPH elastic dynamics, journal of Computer methods
in applied mechanics and engineering (2001)
"""
# Dim of Rings
dim = 2

# Fill up with zeros to 3D
fillUp = False

# ring properties: inner and outer radius
r_inner = 3.0
r_outer = 4.0

# particle spacing

# particle spacing
delta_p = 0.2
# 0.1   --> 2 * 2.196   = 4.392   particles
# 0.07  --> 2 * 4.488   = 8.976   particles
# 0.05  --> 2 * 8.804   = 17.608  particles
# 0.045 -->             = 21.704  particles
# 0.04  --> 2 * 13.736  = 27.472  particles
# 0.038 -->             = 30.384  particles
# 0.037 -->             = 32.136  particles
# 0.035 -->             = 35.856  particles
# 0.034 -->             = 38.048  particles
# 0.033 -->             = 40.352  particles
# 0.03  --> 2 * 24.420  = 48.840  particles
# 0.028 -->             = 56.024  particles
# 0.027 -->             = 60.384  particles
# 0.0269 ->             = 60.784
# 0.025 -->             = 70.416  particles
# 0.02  --> 2 * 54.988  = 109.976 particles
# 0.019 -->             = 121.792 particles
# 0.017 -->             = 152.184 particles
# 0.015 --> 2 *         = 195.648 particles
# 0.014 -->             = 224.288
# 0.0135                = 241.176
# 0.01345               = 243.024
# 0.0134                = 245.072
# 0.013 -->             = 260.264 particles
# 0.0125 ->             = 281.472 particles
# 0.012 -->             = 305.368 particles
# 0.011 -->             = 363.568 particles
# 0.01  --> 2 * 219.860 = 439.720 particles
# 0.00951               = 486.248
# 0.0095                = 487.224
# 0.009 -->             = 542.928 particles
# 0.008 --> 2 * 343.668 = 687.336 particles
# 0.007 --> 2 * 448.772 = 897.544 particles
# 0.0068                = 951.104
# 0.006725              = 972.296
# 0.00672               = 973.960
# 0.00671               = 976.832
# 0.0067                = 979.832
# 0.0066 -->            = 1.009.760 particles
# 0.0065 ->             = 1.040.984 particles
# 0.006 --> 2 * 610.948 = 1.221.896 particles
# 0.005 --> 2 * 879.624 = 1.759.248 particles
# 0.004756              = 1.944.504
# 0.004755              = 1.945.112
# 0.004754              = 1.945.864
# 0.004752              = 1.947.640
# 0.004751              = 1.948.616
# 0.00475               = 1.949.248
# 0.0047                = 1.991.176
# 0.0045 -->            = 2.171.832 particles
# 0.0042                = 2.493.408
# 0.004 -->             = 2.749.000 particles
# 0.0039                = 2.891.648
# 0.00385               = 2.967.128
# 0.00384               = 2.982.752
# 0.00383               = 2.998.400
# 0.00382               = 3.013.912
# 0.0038                = 3.045.720
# 0.0037                = 3.212.648
# 0.0035 ->             = 3.590.312 particles
# 0.003 -->             = 4.886.560 particles
# 0.00251               = 6.981.016 particles
# 0.0025 ->             = 7.037.216 particles
# 0.00249 >             = 7.093.392 particles
# 0.00248 >             = 7.150.944 particles
# 0.002475              = 7.180.096 particles
# 0.0023                = 8.314.416 particles
# 0.0021 ->             = 9.973.296 particles
# 0.002095              = 10.020.728 particles
# 0.00209 >             = 10.068.952 particles
# 0.00208 >             = 10.166.272 particles
# 0.00206 >             = 10.364.248 particles
# 0.00205 >             = 10.466.120 particles
# 0.002 -->             = 10.995.528 particles
# 0.0018 ->             = 13.574.592 particles
# 0.00175 >             = 14.361.496 particles
# 0.0014 ->             = 22.439.880 particles
# 0.0013 ->             = 26.024.920 particles
# 0.00125 >             = 28.148.816 particles
# 0.00124 >             = 28.604.048 particles
# 0.001239              = 28.650.568 particles
# 0.001238              = 28.697.560 particles
# 0.0012375             = 28.719.952 particles
# 0.001 --> 2 * 21.991.108 particles
# 0.0009 ->             = 54.299.264 particles
# 0.000899              = 54.419.648
# 0.000898              = 54.540.928
# 0.00089 >             = 55.526.112
# 0.0007 ->             = 89.760.104 particles
# 0.00067 >             = 97.977.872 particles
# 0.00066 >             = 100.969.152 particles
# 0.00065 >             = 104.099.968 particles
# 0.0005 ->             = 175.929.296 particles

# shift of the rings from origin on x-axis
# shift = 6    # for delta_p <= 0.1
shift = 5
# shift = 4.85 # for delta_p <= 0.01, 0.005, 0.001
# shift = 4.75  # for delta_p <= 0.05
# shift = 4.2  # for delta p <= 0.002

# projected speed
v_p = 0.059 # should be for testcase with AS
# v_p = 0.03 # testcase without AS
# v_p = 100 # For testing with mfm

density = 1
mass = delta_p**2 * density

# ----------------------------------------------------------------------
# Constants for specific internal energy (same as in generateIC.py)
P = 2.5                     # pressure constant
gamma = 5.0 / 3.0           # adiabatic index
u_const = P / ((gamma - 1.0) * density)   # = 3.75 for density = 1
# ----------------------------------------------------------------------

# create initial distribution through creating a 2d grid with dims 2*r_outer x 2*r_outer
# then delete particles which are not on the ring

# calc number of particles within square
N_length = int(2*r_outer/delta_p)
# number of particles in square
N_square = int(N_length**2)

# coordinates of particles in square
r = np.zeros((N_square, dim))

# 2D meshgrid
a = np.mgrid[0:N_length, 0:N_length]

# create square
for i in range(dim):
    k = 0
    for j in range(N_square):
        if j % N_length == 0 and j > 0:
            k += 1
            k = k % N_length
        # print(i, k, j)

        r[j, i] = (a[i, k, j % N_length]-(N_length-1)/2)*delta_p

# count particles in one ring
N = 0
arr = np.zeros(N_square) 
for i in range(N_square):
    radius = np.sqrt(r[i, 0]**2 + r[i, 1]**2)
    if r_outer >= radius >= r_inner:
        N += 1
        arr[i] = 1

# construct two rings with N particles which then are shifted along the x-axis
if fillUp:
    r_ring = np.zeros((N, dim+1))  # first ring
    r_ring2 = np.zeros((N, dim+1))  # second ring
    v = np.zeros((N, dim+1))
    v2 = np.zeros((N, dim+1))
else:
    r_ring = np.zeros((N, dim))  # first ring
    r_ring2 = np.zeros((N, dim))  # second ring
    v = np.zeros((N, dim))
    v2 = np.zeros((N, dim))

m = np.ones(2*N) * mass           # masses (same for all particles)
rho = np.ones(2*N) * density      # (kept only for internal use, not written)
materialId = np.zeros(2*N, dtype=np.int8)
u = np.ones(2*N) * u_const        # specific internal energy (constant)

# create ring 1
counter = 0
for i in range(N_square):
    if arr[i] == 1:
        r_ring[counter, 0] = r[i, 0] - shift
        r_ring[counter, 1] = r[i, 1]
        if fillUp:
            r_ring[counter, 2] = 0

        v[counter, 0] = v_p
        counter += 1

# create ring 2
counter = 0
for i in range(N_square):
    if arr[i] == 1:
        r_ring2[counter, 0] = r[i, 0] + shift
        r_ring2[counter, 1] = r[i, 1]
        if fillUp:
            r_ring2[counter, 2] = 0

        v2[counter, 0] = -v_p
        counter += 1

# put two rings in one array
r_final = np.concatenate((r_ring, r_ring2))
v_final = np.concatenate((v, v2))

if fillUp:
    h5f = h5py.File("rings_deltap{}-3D.h5".format(delta_p), "w")
    print("Saving to rings_deltap{}-3D.h5...".format(delta_p))
else:
    h5f = h5py.File("rings_deltap{}-2D.h5".format(delta_p), "w")
    print("Saving to rings_deltap{}-2D.h5...".format(delta_p))

# Write datasets in the style of generateIC.py
h5f.create_dataset("x", data=r_final)
h5f.create_dataset("v", data=v_final)
h5f.create_dataset("m", data=m)
h5f.create_dataset("u", data=u)               # specific internal energy
h5f.create_dataset("materialId", data=materialId)
h5f.create_dataset("time", data=0.0*m)             # time (scalar, set to zero)

h5f.close()
print("Number of particles: ", 2*N)
print("Finished")
