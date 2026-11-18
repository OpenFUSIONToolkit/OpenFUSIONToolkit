!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file mhd_utils.F90
!
!> Plasma parameters and evaluation routines
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
MODULE mhd_utils
USE oft_local
IMPLICIT NONE
!---
REAL(r8), PARAMETER :: mu0 = pi*(4.E-7_r8) !< \f$ \mu_0 \f$ - Permeability of free-space
REAL(r8), PARAMETER :: ep0 = 8.854187817E-12_r8 !< \f$ \epsilon_0 \f$ - Permitivity of free-space
REAL(r8), PARAMETER :: proton_mass = 1.672621777E-27_r8 !< \f$ m_p \f$ - Proton mass [Kg]
REAL(r8), PARAMETER :: elec_mass = 9.10938291E-31_r8 !< \f$ m_e \f$ - Electron mass [Kg]
REAL(r8), PARAMETER :: elec_charge = 1.602176565E-19_r8 !< \f$ q_e \f$ - Electron charge [C]
!---
CONTAINS
!------------------------------------------------------------------------------
!> Evalute Braginskii transport coefficients for a combined temperature eq
!!
!! @note The highest physical conduction coefficient is chosen for each
!! direction (usually parallel->electron, perpendicular->ion)
!------------------------------------------------------------------------------
PURE FUNCTION brag_comb_transport(n,T,mu,bmag,lambda,me_factor) RESULT(chi)
REAL(r8), INTENT(in) :: n !< Density [m^-3]
REAL(r8), INTENT(in) :: T !< Temperature [eV]
REAL(r8), INTENT(in) :: mu !< Reduced ion mass
REAL(r8), INTENT(in) :: bmag !< Magnetic field strength [T]
REAL(r8), OPTIONAL, INTENT(in) :: lambda !< Coulomb logarithm {11} (optional)
REAL(r8), OPTIONAL, INTENT(in) :: me_factor !< Electron mass enhancement {1} (optional)
REAL(r8) :: chi(2) !< \f$ \left( \chi_{\parallel,e}, min(\chi_{\perp,i},\chi_{\parallel,i}) \right) \f$
REAL(r8) :: chitmp(2)
chi=brag_ion_transport(n,T,mu,bmag,lambda)
chitmp=brag_elec_transport(n,T,bmag,lambda,me_factor)
chi(1)=chitmp(1) ! Use electron parallel
chi(2)=MAX(chi(2),chitmp(2)) ! Use electron perp if magnetization is low enough
END FUNCTION brag_comb_transport
!------------------------------------------------------------------------------
!> Evalute Braginskii transport coefficients for the ion fluid
!!
!! @note The highest physical conduction coefficient is chosen for each direction
!! (usually parallel->electron, perpendicular->ion)
!------------------------------------------------------------------------------
PURE FUNCTION brag_ion_transport(n,T,mu,bmag,lambda) RESULT(chi)
REAL(r8), INTENT(in) :: n !< Density [m^-3]
REAL(r8), INTENT(in) :: T !< Temperature [eV]
REAL(r8), INTENT(in) :: mu !< Reduced ion mass
REAL(r8), INTENT(in) :: bmag !< Magnetic field strength [T]
REAL(r8), OPTIONAL, INTENT(in) :: lambda !< Coulomb logarithm (optional: default = 11)
REAL(r8) :: chi(2) !< \f$ \left( \chi_{\parallel,i}, min(\chi_{\perp,i},\chi_{\parallel,i}) \right) \f$
REAL(r8) :: lam,frac,mag_ion
lam=11.d0
IF(PRESENT(lambda))lam=lambda
!---Compute magnetization parameter
mag_ion=elec_charge*bmag*2.09d13*(T**1.5d0)/(n*lam*proton_mass*SQRT(mu))
!---Compute ion conduction coefficients
chi(1)=(3.906d0*2.09d13)*elec_charge*(T**2.5d0)/(SQRT(mu)*proton_mass*lam) ! Parallel
frac=(2.d0*mag_ion**2 + 2.645d0)/(mag_ion**4 + 2.7d0*mag_ion**2 + 0.677d0)
chi(2)=2.09d13*elec_charge*(T**2.5d0)*frac/(SQRT(mu)*proton_mass*lam)
chi(2)=MIN(chi(1),chi(2)) ! Ensure limit is parallel
END FUNCTION brag_ion_transport
!------------------------------------------------------------------------------
!> Evalute Braginskii transport coefficients for the electron fluid
!
! in our units we have 6*sqrt(2*pi)*pi*ep0^2*sqrt(m_ion)*k_boltz^(5/2)/q^5
! (with k_boltz = elec_charge) this = 2.09e13
!------------------------------------------------------------------------------
PURE FUNCTION brag_elec_transport(n,T,bmag,lambda,me_factor) RESULT(chi)
REAL(r8), INTENT(in) :: n !< Density [m^-3]
REAL(r8), INTENT(in) :: T !< Temperature [eV]
REAL(r8), INTENT(in) :: bmag !< Magnetic field strength [T]
REAL(r8), OPTIONAL, INTENT(in) :: lambda !< Coulomb logarithm (optional: default = 11)
REAL(r8), OPTIONAL, INTENT(in) :: me_factor !< Electron mass enhancement {1} (optional)
REAL(r8) :: chi(2) !< \f$ \left( \chi_{\parallel,e}, min(\chi_{\perp,e},\chi_{\parallel,e}) \right) \f$
REAL(r8) :: lam,me_fac,frac,mag_elec
lam=11.d0
IF(PRESENT(lambda))lam=lambda
me_fac=1.d0
IF(PRESENT(me_factor))me_fac=me_factor
!---Compute magnetization parameter
mag_elec=elec_charge*bmag*3.44d11*(T**1.5d0)/(n*lam*elec_mass*SQRT(me_fac))
!---Compute electron conduction coefficients
chi(1)=(3.1616d0*3.44d11)*elec_charge*(T**2.5d0)/(SQRT(me_fac)*elec_mass*lam) ! Parallel
frac=(4.664d0*mag_elec**2 + 11.92d0)/(mag_elec**4 + 14.79d0*mag_elec**2 + 3.7703d0)
chi(2)=3.44d11*elec_charge*(T**2.5d0)*frac/(SQRT(me_fac)*elec_mass*lam) ! Perpendicular
chi(2)=MIN(chi(2),chi(1)) ! Ensure limit is parallel
END FUNCTION brag_elec_transport
!------------------------------------------------------------------------------
!> Evalute electron-ion thermalization rate
!!
!! Definition from page 34 of NRL plasma formulary
!------------------------------------------------------------------------------
PURE FUNCTION elec_ion_therm_rate(Te,n,mu) RESULT(tau_th)
REAL(r8), INTENT(in) :: Te !< Electron temperature [eV]
REAL(r8), INTENT(in) :: n !< Density [m^-3]
REAL(r8), INTENT(in) :: mu !< Reduced ion mass
REAL(r8) :: tau_th !< \f$ \tau = 1 / (3.2E-15 * Z^2 \lambda / (\mu * T^{3/2})) \f$
tau_th = 3.125d14*mu*(Te**1.5d0)/(log_lambda(Te,n)*n)
END FUNCTION elec_ion_therm_rate
!------------------------------------------------------------------------------
!> Evalute ln( Lambda )
!!
!! Definition for electron-ion collisions in regime \f$ T_i m_e / m_i < 10 Z^2 < T_e \f$,
!! see page 34 of NRL plasma formulary
!------------------------------------------------------------------------------
PURE FUNCTION log_lambda(T,n) RESULT(lambda)
REAL(r8), INTENT(in) :: T !< Temperature [eV]
REAL(r8), INTENT(in) :: n !< Density [m^-3]
REAL(r8) :: lambda !< \f$ \lambda = 24 - ln( \sqrt{n/10^6} / T ) \f$
lambda=24.d0-LOG(SQRT(n/1.d6)/T)
END FUNCTION log_lambda
!------------------------------------------------------------------------------
!> Evalute simple ion viscosity
!!
!! \f$ \tau_i = 2.09 \times 10^7 T_i^{3/2} \mu^{1/2} / (n/10^6) \lambda \f$
!------------------------------------------------------------------------------
PURE FUNCTION ion_visc(n,T,mu,lam) RESULT(nu)
REAL(r8), INTENT(in) :: n !< Density [m^-3]
REAL(r8), INTENT(in) :: T !< Temperature [eV]
REAL(r8), OPTIONAL, INTENT(in) :: mu !< Reduced ion mass (optional: default = 1.0)
REAL(r8), OPTIONAL, INTENT(in) :: lam !< Coulomb logarithm (optional: default = 11)
REAL(r8) :: nu !< \f$ 0.96 n k T_i \tau_i \f$
REAL(r8) :: tau,lambda,mi
lambda=11.d0
IF(PRESENT(lam))lambda=lam
mi=1.d0
IF(PRESENT(mu))mi=mu
tau = 2.09d7*(T**1.5d0)*SQRT(mi)/(n*lambda/1.d6)
nu = 0.96d0*n*elec_charge*T*tau
END FUNCTION ion_visc
!------------------------------------------------------------------------------
!> Evalute Spitzer resistivity
!------------------------------------------------------------------------------
PURE FUNCTION eta_spitzer(T,lam) RESULT(eta)
REAL(r8), INTENT(in) :: T !< Temperature [eV]
REAL(r8), OPTIONAL, INTENT(in) :: lam !< !< Coulomb logarithm (optional: default = 11)
REAL(r8) :: eta !< \f$ 5.253e-05 \frac{\lambda}{T^{3/2} \mu_0} \f$
REAL(r8) :: lambda
lambda=11.d0
IF(PRESENT(lam))lambda=lam
eta = 5.253d-5*lambda/((T**1.5d0)*mu0) ! eta_parallel = eta_perp/1.96
END FUNCTION eta_spitzer
!------------------------------------------------------------------------------
!> Evalute Chodura resistivity
!------------------------------------------------------------------------------
PURE FUNCTION eta_chodura(j,T,n,m,m_factor) RESULT(eta)
REAL(r8), INTENT(in) :: j !<  Current density [A/m^2]
REAL(r8), INTENT(in) :: T !< Temperature [eV]
REAL(r8), INTENT(in) :: n !< Number density [m^-3]
REAL(r8), INTENT(in) :: m !< Ion mass [kg]
REAL(r8), OPTIONAL, INTENT(in) :: m_factor !< Electron mass enhancement (optional)
REAL(r8) :: eta !< \f$ \frac{C m_e \omega_{p,i}}{n e^2} ( 1 - e^{-v_{d,e} / f v_{s,i}} ) \f$
REAL(r8) :: lam,mf
REAL(r8), PARAMETER :: C = 0.1d0
REAL(r8), PARAMETER :: f = 3.d0
mf=1.d0
IF(PRESENT(m_factor))mf=m_factor
!---
eta = C*elec_mass*mf*plasma_freq(n,m)*(1.d0 &
- EXP(-electron_drift_speed(j,n)/(f*sound_speed(T,m))))/(n*(elec_charge**2))
END FUNCTION eta_chodura
!------------------------------------------------------------------------------
!> Evalute plasma frequency for singly charged species
!------------------------------------------------------------------------------
PURE FUNCTION plasma_freq(n,m) RESULT(wp)
REAL(r8), INTENT(in) :: n !< Number density [m^-3]
REAL(r8), INTENT(in) :: m !< Particle mass [kg]
REAL(r8) :: wp !< \f$ \left( 4 \pi n e^2 / m \right)^{1/2} \f$
REAL(r8) :: lam
wp = SQRT(n*(elec_charge**2)/(m*ep0))
END FUNCTION plasma_freq
!------------------------------------------------------------------------------
!> Evalute sound speed
!------------------------------------------------------------------------------
PURE FUNCTION sound_speed(T,m) RESULT(cs)
REAL(r8), INTENT(in) :: T !< Fluid temperature [eV]
REAL(r8), INTENT(in) :: m !< Particle mass [kg]
REAL(r8) :: cs !< \f$ \left( 3 k T / m \right)^{1/2} \f$
REAL(r8), PARAMETER :: gamma = 3.d0
cs = SQRT(gamma*elec_charge*T/m)
END FUNCTION sound_speed
!------------------------------------------------------------------------------
!> Evalute drift velocity for Electrons
!------------------------------------------------------------------------------
PURE FUNCTION electron_drift_speed(j,n) RESULT(vd)
REAL(r8), INTENT(in) :: j !< Current density [A/m^2]
REAL(r8), INTENT(in) :: n !< Number density [m^-3]
REAL(r8) :: vd !< \f$ j / ne \f$
vd = j/(n*elec_charge)
END FUNCTION electron_drift_speed
END MODULE mhd_utils
