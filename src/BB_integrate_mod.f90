!!!
! Elspeth KH Lee - Nov 2026 : Shared band-integrated Planck function for the corr-k lw schemes
!                - Nov 2026 : Optimisation - series exponentials exp(-j x) as powers of exp(-x)
!                             (one exponential per band edge), reciprocals 1/j tabulated
!
! Integrates the blackbody intensity between two wavenumbers using the sum convergence method
! taken from: spectralcalc.com/blackbody/inband_radiance.html (Widger & Woodall 1976)
!!!

module BB_integrate_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: hp = 6.62607015e-34_dp
  real(dp), parameter :: kb = 1.380649e-23_dp
  real(dp), parameter :: c_s = 2.99792458e8_dp
  real(dp), parameter :: c1 = (hp * c_s) / kb
  real(dp), parameter :: c2 = c_s**2
  real(dp), parameter :: n2 = 2.0_dp * hp * c2

  !! Floor value for the band intensity - the lw kernels inherited from the semi-grey code
  !! (e.g. AA_E, VIM_E) take logarithms of Planck ratios, which require be > 0
  real(dp), parameter :: be_floor = 1.0e-99_dp

  !! Maximum number of series terms, the exponent cap of each term (exp(-min(j x, xcap))),
  !! and the reciprocals 1/j
  integer, parameter :: nterm = 151
  real(dp), parameter :: xcap = 300.0_dp
  integer, private :: jj
  real(dp), dimension(nterm), parameter :: rj = (/(1.0_dp/real(jj,dp), jj = 1, nterm)/)

  public :: BB_integrate

contains

  subroutine BB_integrate(n_b, Te, wn_e, be)
    implicit none

    integer, intent(in) :: n_b
    real(dp), intent(in) :: Te
    real(dp), dimension(n_b+1), intent(in) :: wn_e

    real(dp), dimension(n_b), intent(out) :: be

    integer :: ww, j, intitera
    real(dp), dimension(n_b+1) :: iB
    real(dp) :: x, x2, x3, itera, summ, dn, e1, ej, ecap

    !! Zero temperature (e.g. Tint = 0) has no emission
    if (Te < 1e-6_dp) then
      be(:) = 0.0_dp
      return
    end if

    ecap = exp(-xcap)

    do ww = 1, n_b+1

      x = c1 * 100.0_dp * wn_e(ww)/ Te
      x2 = x**2
      x3 = x**3

      itera = 2.0_dp + 20.0_dp/x
      if (itera > 150) then
        itera = 150
      end if
      intitera = int(itera)

      !! exp(-min(j x, xcap)) = exp(-x)**j until j x passes the cap
      e1 = exp(-min(x, xcap))
      ej = 1.0_dp
      summ = 0.0_dp
      do j = 1, intitera + 1
        dn = rj(j)
        if (real(j,dp)*x < xcap) then
          ej = ej*e1
        else
          ej = ecap
        end if
        summ = summ + ej*(x3 + (3.0_dp * x2 + 6.0_dp*(x+dn)*dn)*dn)*dn
      end do

      iB(ww) = n2 * (Te/c1)**(4) * summ
    end do

    do ww = 1, n_b
      be(ww) = max(iB(ww+1) - iB(ww), be_floor)
    end do

  end subroutine BB_integrate

end module BB_integrate_mod
