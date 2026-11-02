!!!
! Elspeth KH Lee - Jun 2024 : Overhauled
!                - Nov 2026 : Rewritten as a true layer-by-layer adding method
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!
! sw: Adding method for shortwave radiation, used in various Earth GCMs (e.g. Mendonça et al. 2015, RRTMG)
!     Each layer has its own delta-Eddington reflection and transmission (Meador & Weaver 1980),
!     then the layers are added together with one upward and one downward pass
!     Pros: Very fast, no matrix inversions, multiple scattering, works for clouds and hazes
!     Cons: Two-stream (Eddington) accuracy within each layer
!!!

module sw_adding_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Maximum exponent value
  real(dp), parameter :: e_max = 100.0_dp

  !! How close to the lambda*mu0 = 1 resonance before mu0 is nudged
  real(dp), parameter :: res_eps = 1.0e-6_dp

  private :: sw_adding_method, layer_RT
  public :: sw_adding

contains

  subroutine sw_adding(nlay, nlev, nb, ng, gw, tau_e, mu_z, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev, nb, ng
    real(dp), dimension(ng), intent(in) :: gw
    real(dp), dimension(ng,nb,nlev), intent(in) :: tau_e
    real(dp), dimension(ng,nb,nlay), intent(in) :: ssa, gg
    real(dp), dimension(nlev), intent(in) :: mu_z
    real(dp), dimension(nb), intent(in) :: Finc

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down, sw_net
    real(dp), intent(out) :: asr

    !! Work variables
    integer :: b, g
    real(dp), dimension(nlev) :: sw_up_g, sw_down_g

    sw_up(:) = 0.0_dp
    sw_down(:) = 0.0_dp

    !! Shortwave flux calculation - sum over each band and g-ordinance
    if (mu_z(nlev) > 0.0_dp) then
      do b = 1, nb
        if (Finc(b) <= 0.0_dp) cycle
        do g = 1, ng
          call sw_adding_method(nlay, nlev, Finc(b), mu_z(:), tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), &
            & sw_down_g(:), sw_up_g(:))
          sw_up(:) = sw_up(:) + sw_up_g(:) * gw(g)
          sw_down(:) = sw_down(:) + sw_down_g(:) * gw(g)
        end do
      end do
    end if

    !! Net sw flux
    sw_net(:) = sw_up(:) - sw_down(:)

    !! Absorbed Stellar Radiation (ASR)
    asr = sw_down(1) - sw_up(1)

  end subroutine sw_adding

  subroutine sw_adding_method(nlay, nlev, F0_in, mu_in, tau_in, w_in, g_in, flx_down, flx_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0_in
    real(dp), dimension(nlev), intent(in) :: tau_in, mu_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    !! Work variables
    integer :: k
    real(dp) :: f, den
    real(dp), dimension(nlay) :: dtau, w0, hg, mu_l, T_dir
    real(dp), dimension(nlay) :: R_b, T_b, R_d, T_d
    real(dp), dimension(nlev) :: T_beam, cum_trans
    real(dp), dimension(nlev) :: Ru_d, Ru_b, Rd_d, D_a
    real(dp) :: D

    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w_in(:) <= 1.0e-12_dp)) then

      ! Check mu_z for spherical corrections
      if (mu_in(nlev) == mu_in(1)) then
        ! No zenith correction, use regular method
        flx_down(:) = F0_in * mu_in(nlev) * exp(-tau_in(:)/mu_in(nlev))
      else
        ! Zenith angle correction, use cumulative transmission function
        cum_trans(1) = tau_in(1)/mu_in(1)
        do k = 1, nlev-1
          cum_trans(k+1) = cum_trans(k) + (tau_in(k+1) - tau_in(k))/mu_in(k+1)
        end do
        flx_down(:) = F0_in * mu_in(nlev) * exp(-cum_trans(:))
      end if

      flx_up(:) = 0.0_dp

      return

    end if

    !! Delta-Eddington scaling (f = g^2)
    do k = 1, nlay
      f = g_in(k)**2
      dtau(k) = (1.0_dp - w_in(k)*f) * (tau_in(k+1) - tau_in(k))
      w0(k) = w_in(k) * (1.0_dp - f)/(1.0_dp - w_in(k)*f)
      hg(k) = g_in(k)/(1.0_dp + g_in(k))
    end do

    !! Zenith angle in each layer and direct beam transmission to each level
    if (mu_in(nlev) == mu_in(1)) then
      mu_l(:) = mu_in(nlev)
    else
      mu_l(:) = mu_in(2:nlev)
    end if
    cum_trans(1) = 0.0_dp
    do k = 1, nlay
      T_dir(k) = exp(-min(dtau(k)/mu_l(k), e_max))
      cum_trans(k+1) = cum_trans(k) + dtau(k)/mu_l(k)
    end do
    T_beam(:) = exp(-min(cum_trans(:), e_max))

    !! Reflection and transmission of each layer
    do k = 1, nlay
      call layer_RT(dtau(k), w0(k), hg(k), mu_l(k), T_dir(k), R_b(k), T_b(k), R_d(k), T_d(k))
    end do

    !! Upward pass - reflection of everything below each level for the beam (Ru_b) and diffuse (Ru_d)
    !! Gas giant lower boundary - nothing is reflected back up
    Ru_b(nlev) = 0.0_dp
    Ru_d(nlev) = 0.0_dp
    do k = nlay, 1, -1
      den = 1.0_dp/(1.0_dp - R_d(k)*Ru_d(k+1))
      Ru_d(k) = R_d(k) + T_d(k)**2 * Ru_d(k+1) * den
      Ru_b(k) = R_b(k) + T_d(k) * (T_dir(k)*Ru_b(k+1) + T_b(k)*Ru_d(k+1)) * den
    end do

    !! Downward pass - reflection of everything above each level (Rd_d) and the diffuse
    !! downward flux from the layers above (D_a)
    Rd_d(1) = 0.0_dp
    D_a(1) = 0.0_dp
    do k = 1, nlay
      den = 1.0_dp/(1.0_dp - Rd_d(k)*R_d(k))
      ! Diffuse down at the top of layer k, including reflections with the layers above
      D = (D_a(k) + Rd_d(k)*T_beam(k)*R_b(k)) * den
      D_a(k+1) = T_beam(k)*T_b(k) + T_d(k)*D
      Rd_d(k+1) = R_d(k) + T_d(k)**2 * Rd_d(k) * den
    end do

    !! Combine the layers above and below each level
    do k = 1, nlev
      D = (D_a(k) + Rd_d(k)*T_beam(k)*Ru_b(k))/(1.0_dp - Rd_d(k)*Ru_d(k))
      flx_down(k) = T_beam(k) + D
      flx_up(k) = T_beam(k)*Ru_b(k) + D*Ru_d(k)
    end do

    !! Scale with the incident flux
    flx_down(:) = flx_down(:) * mu_in(nlev) * F0_in
    flx_up(:) = flx_up(:) * mu_in(nlev) * F0_in

  end subroutine sw_adding_method

  !! Delta-Eddington reflection and transmission of a single layer (as RRTMG reftra)
  !! R_b, T_b - beam reflection and diffuse transmission, R_d, T_d - diffuse reflection and transmission
  !! Only decaying exponentials are used, so thick layers are stable
  subroutine layer_RT(dtau, w0, g0, mu0, T_dir, R_b, T_b, R_d, T_d)
    implicit none

    real(dp), intent(in) :: dtau, w0, g0, mu0, T_dir
    real(dp), intent(out) :: R_b, T_b, R_d, T_d

    real(dp) :: g1, g2, g3, g4, a1, a2, lam, mu, lm, kmu
    real(dp) :: e1, e2, den, Tm, r1, r2, r3, r4, r5, t1, t2, t3

    !! Eddington coefficients
    g1 = (7.0_dp - w0*(4.0_dp + 3.0_dp*g0))/4.0_dp
    g2 = -(1.0_dp - w0*(4.0_dp - 3.0_dp*g0))/4.0_dp
    !! lambda (factored form avoids cancellation as w0 -> 1)
    lam = sqrt(3.0_dp*(1.0_dp - w0)*(1.0_dp - w0*g0))

    !! Diffuse reflection and transmission
    e1 = exp(-min(lam*dtau, e_max))
    e2 = e1**2
    if (lam*dtau < 1.0e-8_dp) then
      ! Conservative scattering limit
      R_d = g1*dtau/(1.0_dp + g1*dtau)
      T_d = 1.0_dp/(1.0_dp + g1*dtau)
    else
      den = (lam + g1) + (lam - g1)*e2
      R_d = g2*(1.0_dp - e2)/den
      T_d = 2.0_dp*lam*e1/den
    end if

    !! Direct beam - nudge mu0 off the lambda*mu0 = 1 resonance (0/0 in the formula)
    mu = mu0
    kmu = lam*mu
    if (abs(1.0_dp - kmu) < res_eps) then
      mu = mu0 + 2.0_dp*res_eps/max(lam, 1.0e-30_dp)
      kmu = lam*mu
    end if
    g3 = (2.0_dp - 3.0_dp*g0*mu)/4.0_dp
    g4 = 1.0_dp - g3

    !! Conservative scattering limit (Meador & Weaver 1980): the general formulas are 0/0 at lam = 0
    if (lam*dtau < 1.0e-8_dp) then
      R_b = (g1*dtau + (g3 - g1*mu0)*(1.0_dp - T_dir))/(1.0_dp + g1*dtau)
      T_b = 1.0_dp - R_b - T_dir
      R_b = max(R_b, 0.0_dp)
      T_b = max(T_b, 0.0_dp)
      return
    end if
    a1 = g1*g4 + g2*g3
    a2 = g1*g3 + g2*g4
    lm = 1.0_dp - kmu**2

    r1 = (1.0_dp - kmu)*(a2 + lam*g3)
    r2 = (1.0_dp + kmu)*(a2 - lam*g3)
    r3 = 2.0_dp*lam*(g3 - a2*mu)
    r4 = lm*(lam + g1)
    r5 = lm*(lam - g1)
    t1 = (1.0_dp + kmu)*(a1 + lam*g4)
    t2 = (1.0_dp - kmu)*(a1 - lam*g4)
    t3 = 2.0_dp*lam*(g4 + a1*mu)

    ! Beam transmission with the same (nudged) mu so the resonant terms cancel
    Tm = T_dir
    if (mu /= mu0) Tm = exp(-min(dtau/mu, e_max))

    R_b = w0*(r1 - r2*e2 - r3*e1*Tm)/(r4 + r5*e2)
    T_b = -w0*((t1 - t2*e2)*Tm - t3*e1)/(r4 + r5*e2)

    R_b = max(R_b, 0.0_dp)
    T_b = max(T_b, 0.0_dp)

  end subroutine layer_RT

end module sw_adding_mod
