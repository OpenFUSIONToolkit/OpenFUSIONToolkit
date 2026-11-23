import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
#
m = 2.
q = 1.9 #1.2 #1.05
Bz = 0.8
wl = 2.5E-3
L = 0.8
ra = 0.1
Te = 50.0
n = 1000
mu0 = np.pi*4.E-7
rho0 = 2.E19*1.672621777E-27
B0 = 0.1*Bz
r = np.linspace(0.,3*ra,n)
#J = np.ones(r.shape)
# J = 0.5*(np.tanh((ra-r)/wl)+1.)
J = 1.0-np.power(r/ra,6); J[J<0.0] = 0.0
Bp = np.zeros(J.shape)
Bp[1:] = integrate.cumulative_trapezoid(r*J,r)/r[1:]
rho_p = 1.
rho_v = 1.E-1
# rho = 0.5*(rho_p-rho_v)*(np.tanh((ra-r)/wl)+1.) + rho_v
rho = np.zeros(J.shape)
rho[:-1] = integrate.cumulative_trapezoid(J*Bp,r)/(1.602E-19*2.0*Te)/mu0
rho = -(rho-rho[-2]); rho[-1]=0.0; rho += rho[0]*rho_v/2.0
sfac = (2.*np.pi*r[5]*Bz)/(Bp[5]*L*q)
Bp = Bp*sfac
q_r = (2.*np.pi*r[1:]*Bz)/(Bp[1:]*L)
sfac2 = np.interp(ra,r[1:],q_r)/q
Bp = Bp*sfac2
q_r /= sfac2
print(q_r[0],np.interp(ra,r[1:],q_r),sfac)
tmp = 1./np.power(1. - (m - q),1./(2.*m))
print('R_w    = {0:.4E}'.format(tmp))
V_a = B0/np.sqrt(mu0*rho0)
V_a_max = B0/np.sqrt(mu0*rho0*rho_v)
print('V_a    = {0:.4E}'.format(V_a))
print('T_a    = {0:.4E} {1:.4E} {2:.4E}'.format(L/V_a,1.733E-02/2.0/V_a,1.733E-02/2.0/V_a_max))
eta = 25.0*np.power(6.0/10.0,3.0/2.0)
print('S      = {0:.4E} {1:.4E}'.format(L*V_a/eta,V_a/(5.E5)))
R = L/(2.*np.pi)
rw = 2.0
gamma = np.sqrt(2.*B0*B0*(m-q)*(1.-np.power(1./rw,2*m) - m + q)/(mu0*rho0*R*R*q*q)/(1.-np.power(1./rw,2*m)))
print('gamma  = {0:.4E} {1:.4E}'.format(gamma, gamma*R/V_a))
print('dt_CFL = {0:.4E}'.format(8.E-3/V_a/3.))
B_max = np.sqrt(Bp.max()*Bp.max() + Bz*Bz)
# print(B_max/Bz)
fig, ax = plt.subplots(1,1)
ax2=ax.twinx()
ax.plot(r,J,'r')
ax2.plot(r,rho,'b--')
ax.plot(r[1:],q_r,'g')
ax.axvline(ra-wl/2.,color='k',linestyle='--')
ax.axvline(ra+wl/2.,color='k',linestyle='--')
plt.show()
#
with open('rwm_eq.dat', 'w+') as fid:
    fid.write('{0}\n'.format(n))
    for i in range(n):
        fid.write('{0:3.12E} {1:3.12E} {2:3.12E} {3:3.12E}\n'.format(r[i], Bp[i], Bz, rho[i]))
