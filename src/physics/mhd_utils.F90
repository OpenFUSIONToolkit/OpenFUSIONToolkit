!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file mhd_utils.F90
!
!> Plasma parameters and evaluation routines
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup doxy_oft_physics
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
! FUNCTION brag_comb_transport
!---------------------------------------------------------------------------
!> Evalute Braginskii transport coefficients for a combined temperature eq
!!
!! @note The highest physical conduction coefficient is chosen for each
!! direction (usually parallel->electron, perpendicular->ion).
!!
!! @param[in] n Density [m^-3]
!! @param[in] T Temperature [eV]
!! @param[in] mu Reduced ion mass
!! @param[in] bmag Magnetic field strength [T]
!! @param[in] lambda Coulomb logarithm {11} (optional)
!! @param[in] me_factor Electron mass enhancement {1} (optional)
!! @result \f$ \left( \chi_{\parallel,e}, min(\chi_{\perp,i},\chi_{\parallel,i}) \right) \f$
!---------------------------------------------------------------------------
PURE FUNCTION brag_comb_transport(n,T,mu,bmag,lambda,me_factor) RESULT(chi)
REAL(r8), INTENT(in) :: n,T,mu,bmag
REAL(r8), OPTIONAL, INTENT(in) :: lambda,me_factor
REAL(r8) :: chi(2),chitmp(2)
chi=brag_ion_transport(n,T,mu,bmag,lambda)
chitmp=brag_elec_transport(n,T,bmag,lambda,me_factor)
chi(1)=chitmp(1) ! Use electron parallel
chi(2)=MAX(chi(2),chitmp(2)) ! Use electron perp if magnetization is low enough
END FUNCTION brag_comb_transport
!---------------------------------------------------------------------------
! FUNCTION brag_ion_transport
!---------------------------------------------------------------------------
!> Evalute Braginskii transport coefficients for the ion fluid
!!
!! @note The highest physical conduction coefficient is chosen for each direction
!! (usually parallel->electron, perpendicular->ion).
!!
!! @param[in] n Density [m^-3]
!! @param[in] T Temperature [eV]
!! @param[in] mu Reduced ion mass
!! @param[in] bmag Magnetic field strength [T]
!! @param[in] lambda Coulomb logarithm {11} (optional)
!! @result \f$ \left( \chi_{\parallel,i}, min(\chi_{\perp,i},\chi_{\parallel,i}) \right) \f$
!---------------------------------------------------------------------------
PURE FUNCTION brag_ion_transport(n,T,mu,bmag,lambda) RESULT(chi)
REAL(r8), INTENT(in) :: n,T,mu,bmag
REAL(r8), OPTIONAL, INTENT(in) :: lambda
REAL(r8) :: chi(2),lam,frac,mag_ion
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
!---------------------------------------------------------------------------
! FUNCTION brag_elec_transport
!---------------------------------------------------------------------------
!> Evalute Braginskii transport coefficients for the electron fluid
!!
!! @param[in] n Density [m^-3]
!! @param[in] T Temperature [eV]
!! @param[in] bmag Magnetic field strength [T]
!! @param[in] lambda Coulomb logarithm {11} (optional)
!! @param[in] me_factor Electron mass enhancement {1} (optional)
!! @result \f$ \left( \chi_{\parallel,e}, min(\chi_{\perp,e},\chi_{\parallel,e}) \right) \f$
! in our units we have 6*sqrt(2*pi)*pi*ep0^2*sqrt(m_ion)*k_boltz^(5/2)/q^5
! (with k_boltz = elec_charge) this = 2.09e13
!---------------------------------------------------------------------------
PURE FUNCTION brag_elec_transport(n,T,bmag,lambda,me_factor) RESULT(chi)
REAL(r8), INTENT(in) :: n,T,bmag
REAL(r8), OPTIONAL, INTENT(in) :: lambda,me_factor
REAL(r8) :: chi(2),lam,me_fac,frac,mag_elec
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
!---------------------------------------------------------------------------
! FUNCTION elec_ion_therm_rate
!---------------------------------------------------------------------------
!> Evalute electron-ion thermalization rate
!!
!! Definition from page 34 of NRL plasma formulary
!!
!! @param[in] Te Electron temperature [eV]
!! @param[in] n Density [m^-3]
!! @param[in] mu Reduced ion mass
!! @result \f$ \tau = 1 / (3.2E-15 * Z^2 \lambda / (\mu * T^{3/2})) \f$
!---------------------------------------------------------------------------
PURE FUNCTION elec_ion_therm_rate(Te,n,mu) RESULT(tau_th)
REAL(r8), INTENT(in) :: Te,n,mu
REAL(r8) :: tau_th
tau_th = 3.125d14*mu*(Te**1.5d0)/(log_lambda(Te,n)*n)
END FUNCTION elec_ion_therm_rate
!---------------------------------------------------------------------------
! FUNCTION log_lambda
!---------------------------------------------------------------------------
!> Evalute ln( Lambda )
!!
!! Definition for electron-ion collisions in regime \f$ T_i m_e / m_i < 10 Z^2 < T_e \f$,
!! see page 34 of NRL plasma formulary
!!
!! @param[in] T Temperature [eV]
!! @param[in] n Density [m^-3]
!! @result \f$ \lambda = 24 - ln( \sqrt{n/10^6} / T ) \f$
!---------------------------------------------------------------------------
PURE FUNCTION log_lambda(T,n) RESULT(lambda)
REAL(r8), INTENT(in) :: T
REAL(r8), INTENT(in) :: n
REAL(r8) :: lambda
lambda=24.d0-LOG(SQRT(n/1.d6)/T)
END FUNCTION log_lambda
!---------------------------------------------------------------------------
! FUNCTION ion_visc
!---------------------------------------------------------------------------
!> Evalute simple ion viscosity
!!
!! \f$ \tau_i = 2.09 \times 10^7 T_i^{3/2} \mu^{1/2} / (n/10^6) \lambda \f$
!!
!! @param[in] n Density [m^-3]
!! @param[in] T Temperature [eV]
!! @param[in] mu Ion mass [au] (1)
!! @param[in] lam \lambda (10)
!! @result \f$ 0.96 n k T_i \tau_i \f$
!---------------------------------------------------------------------------
PURE FUNCTION ion_visc(n,T,mu,lam) RESULT(nu)
REAL(r8), INTENT(in) :: n,T
REAL(r8), OPTIONAL, INTENT(in) :: mu,lam
REAL(r8) :: tau,lambda,mi,nu
lambda=10.d0
IF(PRESENT(lam))lambda=lam
mi=1.d0
IF(PRESENT(mu))mi=mu
tau = 2.09d7*(T**1.5d0)*SQRT(mi)/(n*lambda/1.d6)
nu = 0.96d0*n*elec_charge*T*tau
END FUNCTION ion_visc
!---------------------------------------------------------------------------
! FUNCTION eta_spitzer
!---------------------------------------------------------------------------
!> Evalute Spitzer resistivity
!!
!! @param[in] T Temperature [eV]
!! @param[in] lam \lambda (10)
!! @result \f$ 5.253e-05 \frac{\lambda}{T^{3/2} \mu_0} \f$
!---------------------------------------------------------------------------
PURE FUNCTION eta_spitzer(T,lam) RESULT(eta)
REAL(r8), INTENT(in) :: T
REAL(r8), OPTIONAL, INTENT(in) :: lam
REAL(r8) :: eta,lambda
lambda=10.d0
IF(PRESENT(lam))lambda=lam
eta = 5.253d-5*lambda/((T**1.5d0)*mu0)
END FUNCTION eta_spitzer
!---------------------------------------------------------------------------
! FUNCTION eta_chodura
!---------------------------------------------------------------------------
!> Evalute Chodura resistivity
!!
!! @param[in] j Current density [A/m^2]
!! @param[in] T Temperature [eV]
!! @param[in] n Number density [m^-3]
!! @param[in] m Ion mass [kg]
!! @param[in] m_factor Electron mass enhancement (optional)
!! @result \f$ \frac{C m_e \omega_{p,i}}{n e^2} ( 1 - e^{-v_{d,e} / f v_{s,i}} ) \f$
!---------------------------------------------------------------------------
PURE FUNCTION eta_chodura(j,T,n,m,m_factor) RESULT(eta)
REAL(r8), INTENT(in) :: j,T,n,m
REAL(r8), OPTIONAL, INTENT(in) :: m_factor
REAL(r8) :: eta,lam,mf
REAL(r8), PARAMETER :: C = 0.1d0
REAL(r8), PARAMETER :: f = 3.d0
mf=1.d0
IF(PRESENT(m_factor))mf=m_factor
!---
eta = C*elec_mass*mf*plasma_freq(n,m)*(1.d0 &
- EXP(-electron_drift_speed(j,n)/(f*sound_speed(T,m))))/(n*(elec_charge**2))
END FUNCTION eta_chodura
!---------------------------------------------------------------------------
! FUNCTION plasma_freq
!---------------------------------------------------------------------------
!> Evalute plasma frequency for singly charged species
!!
!! @param[in] n Number density [m^-3]
!! @param[in] m Particle mass [kg]
!! @result \f$ \left( 4 \pi n e^2 / m \right)^{1/2} \f$
!---------------------------------------------------------------------------
PURE FUNCTION plasma_freq(n,m) RESULT(wp)
REAL(r8), INTENT(in) :: n,m
REAL(r8) :: wp,lam
wp = SQRT(n*(elec_charge**2)/(m*ep0))
END FUNCTION plasma_freq
!---------------------------------------------------------------------------
! FUNCTION sound_speed
!---------------------------------------------------------------------------
!> Evalute sound speed
!!
!! @param[in] T Fluid temperature [eV]
!! @param[in] m Particle mass [kg]
!! @result \f$ \left( 3 k T / m \right)^{1/2} \f$
!---------------------------------------------------------------------------
PURE FUNCTION sound_speed(T,m) RESULT(cs)
REAL(r8), INTENT(in) :: T,m
REAL(r8) :: cs
REAL(r8), PARAMETER :: gamma = 3.d0
cs = SQRT(gamma*elec_charge*T/m)
END FUNCTION sound_speed
!---------------------------------------------------------------------------
! FUNCTION electron_drift_speed
!---------------------------------------------------------------------------
!> Evalute drift velocity for Electrons
!!
!! @param[in] j Current density [A/m^2]
!! @param[in] n Number density [m^-3]
!! @result \f$ j / ne \f$
!---------------------------------------------------------------------------
PURE FUNCTION electron_drift_speed(j,n) RESULT(vd)
REAL(r8), INTENT(in) :: j,n
REAL(r8) :: vd
vd = j/(n*elec_charge)
END FUNCTION electron_drift_speed
END MODULE mhd_utils
