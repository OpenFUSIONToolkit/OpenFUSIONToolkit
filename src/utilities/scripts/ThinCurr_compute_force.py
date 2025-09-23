#!/usr/bin/env python
#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
import argparse
import numpy as np
from OpenFUSIONToolkit.io import build_XDMF

def compute_force(J,B,area,use_slow,reg_min,reg_max,rcc_torque=None,torque_cen=None):
    if use_slow:
        F = np.zeros((3,))
        T = np.zeros((3,))
        for j in range(ThinCurr_mesh.nc):
            if (ThinCurr_mesh.reg[j] >= reg_min) and (ThinCurr_mesh.reg[j] <= reg_max):
                Ftmp = np.cross(J[j,:],B[j,:])
                F += area[j]*Ftmp
                if rcc_torque is not None:
                    T += area[j]*np.cross(rcc[j,:]-torque_cen,Ftmp)
    else:
        B = B*area[:,None]
        F_cc = np.cross(J[reg_mask,:],B[reg_mask,:],axis=1)
        F_mean = np.mean(np.sqrt(np.sum(np.power(F_cc,2),axis=1)),axis=0)
        F = np.sum(F_cc,axis=0)
        if rcc_torque is not None:
            T = np.sum(np.cross(rcc_torque[reg_mask,:],np.cross(J[reg_mask,:],B[reg_mask,:])),axis=0)
    # Print results
    F_mean = max(F_mean,np.finfo(F_mean).tiny) # Avoid divide by zero
    if rcc_torque is None:
        print("{0:15.6E} {1:15.6E} {2:15.6E} ({3:15.6E})".format(*F,np.sqrt(np.sum(np.power(F,2)))/F_mean))
    else:
        print("{0:15.6E} {1:15.6E} {2:15.6E} ({3:15.6E}) {4:15.6E} {5:15.6E} {6:15.6E}".format(*F,np.sqrt(np.sum(np.power(F,2)))/F_mean,*T))

parser = argparse.ArgumentParser()
parser.description = "Compute force on regions from ThinCurr simulations"
parser.add_argument("--reg_min", default=-1.0, type=float, help='Minimum region ID')
parser.add_argument("--reg_max", default=1000.0, type=float, help='Maximum region ID')
parser.add_argument("--btr0", default=0.0, type=float, help='Vacuum toroidal flux F=B_t,0 * R_0')
parser.add_argument("--torque_cen", default=None, type=float, nargs='+', help='Torque center')
parser.add_argument("--nmax", default=1000, type=int, help='Max time step or eigenvalue index (if negative)')
parser.add_argument("--use_slow", default=False, action='store_true', help='Use slow surface integration')

options = parser.parse_args()
have_torque = False
if options.torque_cen is not None:
    torque_cen = np.asarray(options.torque_cen)
    have_torque = True
else:
    torque_cen = np.r_[0.0,0.0,0.0]

plot_data = build_XDMF()
ThinCurr_mesh = plot_data['ThinCurr']['smesh']
area = np.zeros((ThinCurr_mesh.nc,))
rcc = np.zeros((ThinCurr_mesh.nc,3))
rcc_torque = np.zeros((ThinCurr_mesh.nc,3))
for i in range(ThinCurr_mesh.nc):
    v1 = ThinCurr_mesh.r[ThinCurr_mesh.lc[i,1],:]-ThinCurr_mesh.r[ThinCurr_mesh.lc[i,0],:]
    v2 = ThinCurr_mesh.r[ThinCurr_mesh.lc[i,2],:]-ThinCurr_mesh.r[ThinCurr_mesh.lc[i,0],:]
    area[i] = np.linalg.norm(np.cross(v1,v2))/2.0
    rcc[i,:] = (ThinCurr_mesh.r[ThinCurr_mesh.lc[i,2],:]+ThinCurr_mesh.r[ThinCurr_mesh.lc[i,1],:]+ThinCurr_mesh.r[ThinCurr_mesh.lc[i,0],:])/3.0
    rcc_torque[i,:] = rcc[i,:]-torque_cen
Btor = np.zeros((ThinCurr_mesh.np,3))
for i in range(ThinCurr_mesh.np):
    that = np.r_[-ThinCurr_mesh.r[i,1],ThinCurr_mesh.r[i,0],0.0]; that/=np.linalg.norm(that)
    Btor[i,:] = options.btr0*that/np.sqrt(np.power(ThinCurr_mesh.r[i,0],2)+np.power(ThinCurr_mesh.r[i,1],2))
if not have_torque:
    rcc_torque = None

if rcc_torque is None:
    print("   Fx              Fy              Fz              |F|/F_rms")
else:
    print("   Fx              Fy              Fz              |F|/F_rms         Tx              Ty              Tz")

reg_mask = np.all((ThinCurr_mesh.reg >= options.reg_min, ThinCurr_mesh.reg <= options.reg_max),axis=0)
if options.nmax > 0:
    for i, time in enumerate(ThinCurr_mesh.times):
        if i >= options.nmax:
            break
        J = ThinCurr_mesh.get_field('J',time)
        Bv = ThinCurr_mesh.get_field('B_v',time) + Btor
        B = (Bv[ThinCurr_mesh.lc[:,0],:]+Bv[ThinCurr_mesh.lc[:,1],:]+Bv[ThinCurr_mesh.lc[:,2],:])/3.0
        compute_force(J,B,area,options.use_slow,options.reg_min,options.reg_max,rcc_torque=rcc_torque,torque_cen=torque_cen)
else:
    for i in range(abs(options.nmax)):
        J = ThinCurr_mesh.static_fields['J{0:02d}'.format(i+1)]
        Bv = ThinCurr_mesh.static_fields['B_v{0:02d}'.format(i+1)]
        B = (Bv[ThinCurr_mesh.lc[:,0],:]+Bv[ThinCurr_mesh.lc[:,1],:]+Bv[ThinCurr_mesh.lc[:,2],:])/3.0
        compute_force(J,B,area,options.use_slow,options.reg_min,options.reg_max,rcc_torque=rcc_torque,torque_cen=torque_cen)
                