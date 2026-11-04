!!!
! Elspeth KH Lee - May 2021 : Initial version
!                - Jan 2022 : Working version
!                - Jun 2024 : Optimisations
!                - Oct 2026 : Exo-FMS_2 version, gas giant lower boundary
!                - Nov 2026 : Resonance guard, zcorr energy conservation fix, renamed sw_Toon89 (was sw_Toon)
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!                - Nov 2026 : Optimisations - one beam exponential per level (shared by the direct beam
!                             and the particular solution), reciprocal mu0, contiguous per
!                             g-ordinance inputs
!
! sw: Two-stream method following the "Toon89" method (Toon et al. 1989)
!     Based on the CHIMERA code by Mike Line, but cleaned up slightly
!     Pros: Fast, decently accurate, multiple scattering method
!     Cons: For combined high ssa and g (0.9+,0.9+) can be unstable and inaccurate
!!!

module sw_Toon89_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Single scattering albedo cap - exactly conservative scattering (w = 1, eigenvalue k = 0) is singular.
  !! 1 - 1e-12 absorbs < 1e-8 of the incident flux even in a tau ~ 6000 column
  real(dp), parameter :: w_max = 1.0_dp - 1.0e-12_dp

  !! How close to the lambda*mu0 = 1 resonance before mu0 is nudged
  real(dp), parameter :: res_eps = 1.0e-6_dp

  private :: sw_Toon89_method, dtridgl
  public :: sw_Toon89

contains

  subroutine sw_Toon89(nlay, nlev, nb, ng, gw, tau_e, mu_z, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
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
    real(dp), dimension(nlev) :: tau_g
    real(dp), dimension(nlay) :: ssa_g, gg_g
    real(dp), parameter :: a_surf = 0.0_dp   ! gas giant: no reflecting lower boundary

    sw_up(:) = 0.0_dp
    sw_down(:) = 0.0_dp

    !! Shortwave flux calculation - sum over each band and g-ordinance
    if (mu_z(nlev) > 0.0_dp) then
      do b = 1, nb
        if (Finc(b) <= 0.0_dp) cycle
        do g = 1, ng
          ! Contiguous copies (avoids array temporaries for the strided slices)
          tau_g(:) = tau_e(g,b,:)
          ssa_g(:) = ssa(g,b,:)
          gg_g(:) = gg(g,b,:)
          call sw_Toon89_method(nlay, nlev, Finc(b), mu_z(:), tau_g, ssa_g, gg_g, a_surf, &
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

  end subroutine sw_Toon89

  subroutine sw_Toon89_method(nlay, nlev, F0_in, mu_in, tau_in, w_in, g_in, w_surf_in, flx_down, flx_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0_in, w_surf_in
    real(dp), dimension(nlev), intent(in) :: tau_in, mu_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    !! Work variables
    integer :: k, i, n
    integer :: l, lm2, lm1

    real(dp), dimension(nlev) :: dir, tau, cum_trans, ct
    real(dp), dimension(nlay) :: dtau_in, dtau, mu_zm, rmu, trans_top, trans_lay, F0_l
    logical, dimension(nlay) :: nudged
    real(dp), dimension(nlay) :: w0, hg
    real(dp), dimension(nlay) :: g1, g2, g3, g4
    real(dp), dimension(nlay) :: lam, gam, denom
    real(dp), dimension(nlay) :: Am, Ap, Cpm1, Cmm1, Cp, Cm
    real(dp), dimension(nlay) :: exptrm, Ep, Em, E1, E2, E3, E4
    real(dp), dimension(nlay+nlay) :: Af, Bf, Cf, Df, xk
    real(dp), dimension(nlay) :: xk1, xk2

    real(dp), parameter :: sqrt3 = sqrt(3.0_dp)
    real(dp), parameter :: sqrt3d2 = sqrt3/2.0_dp
    real(dp), parameter :: bsurf = 0.0_dp, btop = 0.0_dp  ! Surface 'emission' boundary fluxes (0 for shortwave)

    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w_in(:) <= 1.0e-12_dp)) then

      ! Check mu_z for spherical corrections
      if (mu_in(nlev) == mu_in(1)) then
        ! No zenith correction, use regular method
        flx_down(:) = F0_in* mu_in(nlev) * exp(-tau_in(:)/mu_in(nlev))
      else
        ! Zenith angle correction, use cumulative transmission function
        cum_trans(1) = tau_in(1)/mu_in(1)
        do k = 1, nlev-1
          cum_trans(k+1) = cum_trans(k) + (tau_in(k+1) - tau_in(k))/mu_in(k+1)
        end do
        do k = 1, nlev
          flx_down(k) = F0_in * mu_in(nlev) * exp(-cum_trans(k))
        end do
      end if

      flx_down(nlev) = flx_down(nlev) * (1.0_dp - w_surf_in) ! The surface flux for surface heating is the amount of flux absorbed by surface
      flx_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo

      return

    end if

    l = nlay + nlay
    lm2 = l - 2
    lm1 = l - 1

    dtau_in(:) = tau_in(2:nlev) - tau_in(1:nlay)

    !! Delta-Eddington scaling
    w0(:) = ((1.0_dp - g_in(:)**2)*w_in(:))/(1.0_dp-w_in(:)*g_in(:)**2)
    w0(:) = min(w0(:), w_max)
    dtau(:) = (1.0_dp-w_in(:)*g_in(:)**2)*dtau_in(:)
    hg(:) = g_in(:)/(1.0_dp + g_in(:))

    tau(1) = 0.0_dp
    do k = 1, nlay
      tau(k+1) = tau(k) + dtau(k)
    end do

    !! Zenith angle in each layer (zcorr uses the lower level, as the cumulative transmission)
    if (mu_in(nlev) == mu_in(1)) then
      mu_zm(:) = mu_in(nlev)
    else
      mu_zm(:) = mu_in(2:nlev)
    end if

    !! Direct beam transmission to each level
    rmu(:) = 1.0_dp/mu_zm(:)
    cum_trans(1) = tau(1)*rmu(1)
    do k = 1, nlay
      cum_trans(k+1) = cum_trans(k) + dtau(k)*rmu(k)
    end do
    ct(:) = exp(-cum_trans(:))
    dir(:) = F0_in * mu_in(nlev) * ct(:)
    trans_top(:) = ct(1:nlay)

    g1(:) = sqrt3d2 * (2.0_dp-w0(:)*(1.0_dp+hg(:)))
    g2(:) = (sqrt3d2*w0(:)) * (1.0_dp-hg(:))
    where (g2(:) == 0.0_dp)
      g2(:) = 1.0e-10_dp
    end where

    lam(:) = sqrt(g1(:)**2 - g2(:)**2)
    gam(:) = (g1(:) - lam(:))/g2(:)

    !! Nudge mu0 off the lambda*mu0 = 1 resonance (0/0 in the particular solution)
    !! the nudged angle is only used in the layer's own particular solution
    do k = 1, nlay
      nudged(k) = abs(1.0_dp - lam(k)*mu_zm(k)) < res_eps
      if (nudged(k)) then
        mu_zm(k) = mu_zm(k) + 2.0_dp*res_eps/lam(k)
        rmu(k) = 1.0_dp/mu_zm(k)
      end if
    end do
    g3(:) = (1.0_dp - sqrt3*hg(:)*mu_zm(:))/2.0_dp
    g4(:) = 1.0_dp - g3(:)

    !! Beam flux in each layer, so the beam scatters exactly the energy it loses (energy conserving with zcorr)
    F0_l(:) = F0_in * mu_in(nlev)*rmu(:)

    denom(:) = F0_l(:) * w0(:)/(lam(:)**2 - rmu(:)**2)
    Am(:) = denom(:) * (g4(:) * (g1(:) + rmu(:)) + g2(:)*g3(:))
    Ap(:) = denom(:) * (g3(:) * (g1(:) - rmu(:)) + g2(:)*g4(:))

    ! Cpm1 and Cmm1 are the C+ and C- terms evaluated at the top of the layer.
    Cpm1(:) = Ap(:) * trans_top(:)
    Cmm1(:) = Am(:) * trans_top(:)
    ! Cp and Cm are the C+ and C- terms evaluated at the bottom of the layer.
    ! (the beam transmission to the next level, except with the layer's nudged mu0)
    trans_lay(:) = ct(2:nlev)
    do k = 1, nlay
      if (nudged(k)) trans_lay(k) = trans_top(k) * exp(-dtau(k)*rmu(k))
    end do
    Cp(:) = Ap(:) * trans_lay(:)
    Cm(:) = Am(:) * trans_lay(:)

    ! Solve for the coefficients of system of equations using boundary conditions
    ! Exponential terms:
    exptrm(:) = min(lam(:)*dtau(:),35.0_dp)
    Ep(:) = exp(exptrm(:))
    Em(:) = 1.0_dp/Ep(:)

    E1(:) = Ep(:) + gam(:)*Em(:)
    E2(:) = Ep(:) - gam(:)*Em(:)
    E3(:) = gam(:)*Ep(:) + Em(:)
    E4(:) = gam(:)*Ep(:) - Em(:)

    Af(1) = 0.0_dp
    Bf(1) = gam(1) + 1.0_dp
    Cf(1) = gam(1) - 1.0_dp
    Df(1) = btop - Cmm1(1)

    n = 0
    do i = 2, lm2, 2
      n = n + 1
      Af(i) = (E1(n)+E3(n))*(gam(n+1)-1.0_dp)
      Bf(i) = (E2(n)+E4(n))*(gam(n+1)-1.0_dp)
      Cf(i) = 2.0_dp*(1.0_dp-gam(n+1)**2)
      Df(i) = (gam(n+1)-1.0_dp)*(Cpm1(n+1) - Cp(n)) + (1.0_dp-gam(n+1))*(Cm(n)-Cmm1(n+1))
    end do

    n = 0
    do i = 3, lm1, 2
      n = n + 1
      Af(i) = 2.0_dp*(1.0_dp-gam(n)**2)
      Bf(i) = (E1(n)-E3(n))*(1.0_dp + gam(n+1))
      Cf(i) = (E1(n)+E3(n))*(gam(n+1)-1.0_dp)
      Df(i) = E3(n)*(Cpm1(n+1) - Cp(n)) + E1(n)*(Cm(n) - Cmm1(n+1))
    end do

    Af(l) = E1(nlay) - w_surf_in*E3(nlay)
    Bf(l) = E2(nlay) - w_surf_in*E4(nlay)
    Cf(l) = 0.0_dp
    Df(l) = bsurf - Cp(nlay) + w_surf_in*Cm(nlay)

    call dtridgl(l, Af, Bf, Cf, Df, xk)

    do n = 1, nlay
      xk1(n) = xk(2*n-1)+xk(2*n)
      xk2(n) = xk(2*n-1)-xk(2*n)
      if (abs(xk2(n)/xk(2*n-1)) < 1e-30_dp) then
        xk2(n) = 0.0_dp
      end if
    end do

    flx_up(1:nlay) = xk1(:)+gam(:)*xk2(:)+Cpm1(:)
    flx_down(1:nlay) = xk1(:)*gam(:)+xk2(:)+Cmm1(:)

    flx_up(nlev) = xk1(nlay)*Ep(nlay)+gam(nlay)*xk2(nlay)*Em(nlay)+Cp(nlay)
    flx_down(nlev) = xk1(nlay)*Ep(nlay)*gam(nlay)+xk2(nlay)*Em(nlay)+Cm(nlay)

    flx_down(:) = flx_down(:) + dir(:)

  end subroutine sw_Toon89_method

  subroutine dtridgl(l, af, bf, cf, df, xk)
    implicit none

    integer, intent(in) :: l
    real(dp), dimension(l), intent(in) :: af, bf, cf, df
    real(dp), dimension(l), intent(out) :: xk

    integer :: i
    real(dp) :: x
    real(dp), dimension(l) :: as, ds

    as(l) = af(l)/bf(l)
    ds(l) = df(l)/bf(l)

    do i = l-1, 1, -1
      x = 1.0_dp/(bf(i) - cf(i)*as(i+1))
      as(i) = af(i)*x
      ds(i) = (df(i) - cf(i)*ds(i+1))*x
    end do

    xk(1) = ds(1)

    do i = 2, l
      xk(i) = ds(i)-as(i)*xk(i-1)
    end do

  end subroutine dtridgl

end module sw_Toon89_mod
