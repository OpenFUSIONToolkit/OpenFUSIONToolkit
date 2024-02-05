import argparse
import numpy as np
import h5py


def compute_force(J,B,area,use_slow,reg_min,reg_max,rcc_torque=None,torque_cen=None):
    if use_slow:
        F = np.zeros((3,))
        T = np.zeros((3,))
        for j in range(nc):
            if (reg[j] >= reg_min) and (reg[j] <= reg_max):
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

with h5py.File('mesh.0001.h5', 'r') as fid:
    lc = np.asarray(fid['/LC_surf'])
    r = np.asarray(fid['/R_surf'])
nc = lc.shape[0]
area = np.zeros((lc.shape[0],))
rcc = np.zeros((nc,3))
rcc_torque = np.zeros((nc,3))
Btor = np.zeros((nc,3))
for i in range(nc):
    v1 = r[lc[i,1],:]-r[lc[i,0],:]
    v2 = r[lc[i,2],:]-r[lc[i,0],:]
    area[i] = np.linalg.norm(np.cross(v1,v2))/2.0
    rcc[i,:] = (r[lc[i,2],:]+r[lc[i,1],:]+r[lc[i,0],:])/3.0
    rcc_torque[i,:] = rcc[i,:]-torque_cen
    that = np.r_[-rcc[i,1],rcc[i,0],0.0]; that/=np.linalg.norm(that)
    Btor[i,:] = options.btr0*that/np.sqrt(np.power(rcc[i,0],2)+np.power(rcc[i,1],2))
if not have_torque:
    rcc_torque = None

with h5py.File('scalar_dump.0001.h5', 'r') as fid:
    reg = np.asarray(fid['/REG_surf0000'])

reg_mask = np.all((reg >= options.reg_min, reg <= options.reg_max),axis=0)
with h5py.File('vector_dump.0001.h5', 'r') as fid:
    if options.nmax > 0:
        for i in range(options.nmax):
            if '/J{0:04d}'.format(i+1) not in fid:
                break
            J = np.asarray(fid['/J{0:04d}'.format(i+1)])
            B = (np.asarray(fid['/B{0:04d}'.format(i+1)])+Btor)
            compute_force(J,B,area,options.use_slow,options.reg_min,options.reg_max,rcc_torque=rcc_torque,torque_cen=torque_cen)
    else:
        for i in range(abs(options.nmax)):
            if '/J_{0:02d}{1:04d}'.format(i+1,0) not in fid:
                break
            if '/B_{0:02d}{1:04d}'.format(i+1,0) not in fid:
                break
            J = np.asarray(fid['/J_{0:02d}{1:04d}'.format(i+1,0)])
            B = (np.asarray(fid['/B_{0:02d}{1:04d}'.format(i+1,0)])+Btor)
            compute_force(J,B,area,options.use_slow,options.reg_min,options.reg_max,rcc_torque=rcc_torque,torque_cen=torque_cen)
                