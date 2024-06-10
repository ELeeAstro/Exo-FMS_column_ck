!!!
! Elspeth KH Lee - May 2021 : Initial version
!                - Jan 2022 : Working version
!                - Jun 2024 : Optimisations
!
! sw: Two-stream method following the "Toon89" method (Toon et al. 1989)
!     Based on the CHIMERA code by Mike Line, but cleaned up slightly
!     Pros: Fast, decently accurate, multiple scattering method
!     Cons: For combined high ssa and g (0.9+,0.9+) can be unstable and inaccurate
!!!

module sw_Toon_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  private :: sw_Toon89, dtridgl
  public :: sw_Toon

contains

  subroutine sw_Toon(nlay, nlev, nb, ng, gw, tau_e, mu_z, Finc, ssa, gg, a_surf, sw_up, sw_down, sw_net, asr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev, nb, ng
    real(dp), dimension(ng) :: gw
    real(dp), dimension(ng,nb,nlev), intent(in) :: tau_e
    real(dp), dimension(ng,nb,nlay), intent(in) :: ssa, gg
    real(dp), dimension(nlev) :: mu_z
    real(dp), dimension(nb) :: Finc, a_surf

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down, sw_net
    real(dp), intent(out) :: asr

    !! Work variables
    integer :: b, g
    real(dp), dimension(nb,nlev) :: sw_down_b, sw_up_b
    real(dp), dimension(ng,nlev) :: sw_down_g, sw_up_g


    !! Shortwave flux calculation
    if (mu_z(nlev) > 0.0_dp) then
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
      do b = 1, nb
        sw_down_b(b,:) = 0.0_dp
        sw_up_b(b,:) = 0.0_dp
        if (Finc(b) <= 0.0_dp) then
          cycle
        end if
        do g = 1, ng
          call sw_Toon89(nlay, nlev, Finc(b), mu_z(:), tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), a_surf(b), sw_down_g(g,:), sw_up_g(g,:))
          sw_down_b(b,:) = sw_down_b(b,:) + sw_down_g(g,:) * gw(g)
          sw_up_b(b,:) = sw_up_b(b,:) + sw_up_g(g,:) * gw(g)
        end do
        sw_down(:) = sw_down(:) + sw_down_b(b,:)
        sw_up(:) = sw_up(:) + sw_up_b(b,:)
      end do
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Net sw flux
    sw_net(:) = sw_up(:) - sw_down(:)

    !! Absorbed Stellar Radiation (ASR)
    asr = sw_down(1) - sw_up(1)

  end subroutine sw_Toon

  subroutine sw_Toon89(nlay, nlev, F0_in, mu_in, tau_in, w_in, g_in, w_surf_in, flx_down, flx_up)
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

    real(dp), dimension(nlev) :: dir, tau, cum_trans
    real(dp), dimension(nlay) :: dtau_in, dtau, mu_zm
    real(dp), dimension(nlay) :: w0, hg
    real(dp), dimension(nlay) :: g1, g2, g3, g4
    real(dp), dimension(nlay) :: lam, gam, denom
    real(dp), dimension(nlay) :: Am, Ap, Cpm1, Cmm1, Cp, Cm
    real(dp), dimension(nlay) :: exptrm, Ep, Em, E1, E2, E3, E4
    real(dp), dimension(nlay+nlay) :: Af, Bf, Cf, Df, xk
    real(dp), dimension(nlay) :: xk1, xk2

    !! Optimisation Variables
    real(dp), dimension(nlay) :: opt1
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

    do k = 1, nlay
      dtau_in(k) = max(tau_in(k+1) - tau_in(k),1e-6_dp)
    end do

    ! Delta eddington scaling
    w0(:) = ((1.0_dp - g_in(:)**2)*w_in(:))/(1.0_dp-w_in(:)*g_in(:)**2)
    dtau(:) = (1.0_dp-w_in(:)*g_in(:)**2)*dtau_in(:)
    hg(:) = g_in(:)/(1.0_dp + g_in(:))

    tau(1) = 0.0_dp
    do k = 1, nlay
      tau(k+1) = tau(k) + dtau(k)
    end do

    if (mu_in(nlev) == mu_in(1)) then
      ! No zenith correction, use regular method
      dir(:) = F0_in * mu_in(nlev) * exp(-tau(:)/mu_in(nlev))
      mu_zm(:) = mu_in(nlev)
    else
      ! Zenith angle correction, use cumulative transmission function
      cum_trans(1) = tau(1)/mu_in(1)
      do k = 1, nlev-1
        cum_trans(k+1) = cum_trans(k) + (tau(k+1) - tau(k))/mu_in(k+1)
      end do
      do k = 1, nlev
        dir(k) = F0_in * mu_in(nlev) * exp(-cum_trans(k))
      end do
      mu_zm(:) = (mu_in(1:nlay) + mu_in(2:nlev))/2.0_dp ! Zenith angle at midpoints
    end if

    g1(:) = sqrt3d2 * (2.0_dp-w0(:)*(1.0_dp+hg(:)))
    g2(:) = (sqrt3d2*w0(:)) * (1.0_dp-hg(:))
    where (g2(:) <= 1.0e-10_dp)
      g2(:) = 1.0e-10_dp
    end where
    g3(:) = (1.0_dp - sqrt3*hg(:)*mu_zm(:))/2.0_dp
    g4(:) = 1.0_dp - g3(:)

    lam(:) = sqrt(g1(:)**2 - g2(:)**2)
    gam(:) = (g1(:) - lam(:))/g2(:)

    denom(:) = lam(:)**2 - 1.0_dp/(mu_zm(:)**2)
    where (denom(:) <= 1.0e-10_dp)
      denom(:) = 1.0e-10_dp
    end where
    Am(:) = F0_in * w0(:) * (g4(:) * (g1(:) + 1.0_dp/mu_zm(:)) + g2(:)*g3(:))/denom(:)
    Ap(:) = F0_in * w0(:) * (g3(:) * (g1(:) - 1.0_dp/mu_zm(:)) + g2(:)*g4(:))/denom(:)

    ! Cpm1 and Cmm1 are the C+ and C- terms evaluated at the top of the layer.
    opt1(:) = exp(-tau(1:nlay)/mu_zm(:))
    Cpm1(:) = Ap(:) * opt1(:)
    Cmm1(:) = Am(:) * opt1(:)
    ! Cp and Cm are the C+ and C- terms evaluated at the bottom of the layer.
    opt1(:) = exp(-tau(2:nlev)/mu_zm(:))
    Cp(:) = Ap(:) * opt1(:)
    Cm(:) = Am(:) * opt1(:)

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
    flx_up(:) = flx_up(:)

  end subroutine sw_Toon89

  subroutine dtridgl(l, af, bf, cf, df, xk)
    implicit none

    integer, intent(in) :: l
    real(dp), dimension(l), intent(in) :: af, bf, cf, df
    real(dp), dimension(l), intent(out) :: xk

    integer :: i
    integer, parameter :: nmax = 301
    real(dp) :: x
    real(dp), dimension(nmax) :: as, ds

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

end module sw_Toon_mod