!!!
! Elspeth KH Lee - Nov 2026 : Initial version
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!
! lw: delta-two- and four-stream combination method (delta-2/4; Fu et al. 1997, J. Atmos. Sci.),
!     thermal emission only.
!     Step 1:  delta-two-stream discrete ordinates solution at mu1 = 1/sqrt(3) (Fu et al. 1997
!              Eqs. 2.9-2.10 with D = 1/mu1), Planck function linear in tau inside each layer,
!              solved by layer adding (R, T and thermal sources per layer, composites from the top
!              and from the bottom with the gas giant boundary I_up = I_down + B(Tint)).
!     Step 2:  the two-stream solution inside each layer (particular + two exponential modes)
!              gives the scattering source function at the four-stream (double-Gauss) angles,
!              S(+/-mu_i) = (1 - w) B + w/2 [P(+/-mu_i, mu1) I+ + P(+/-mu_i, -mu1) I-],
!              and the RTE is integrated analytically along each four-stream angle
!              (linear and exponential source terms), giving the four-stream intensities and fluxes.
!     The scattering is approximated (two-stream) but the intensity integration is four-stream,
!     which Fu et al. found recovers most of the delta-four-stream accuracy at ~1.8x the
!     two-stream cost (vs ~9x for the full four-stream solution).
!
! Notes:
!   - delta-Eddington scaling (f = g**2) with the two-term phase function P = 1 + 3 g' mu mu':
!     the source function is built from two-stream intensities, so only chi_0 and chi_1 are
!     meaningful. (Tested against delta-M f = g**4 with chi_0..chi_3 retained: the two-term,
!     f = g**2 version was the more accurate against a 48-stream reference.)
!   - With a two-term phase function sum_j [P(mu_i,mu_j) + P(mu_i,-mu_j)]/2 = 1 for any mu1, so
!     energy is conserved exactly; mu1 = 1/1.66 (Elsasser) can be used instead of 1/sqrt(3).
!   - Exact for w = 0 (four-stream linear-in-tau formal solution) and I = B exactly in an
!     isothermal medium.
!   - w is capped at w_max (k -> 0 conservative limit).
!
! Pros: four-stream intensity integration at near two-stream cost, no iterations or matrix inversions
! Cons: scattering still at two-stream accuracy - lw_DDA is more accurate for strongly scattering clouds
!!!

module lw_Fu97_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4
  use BB_integrate_mod, only : BB_integrate
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Two-stream angle (Gauss point 1/sqrt(3); 1/1.66 gives the Elsasser diffusivity version)
  real(dp), parameter :: mu1 = 0.5773502691896258_dp

  !! Four-stream double-Gauss nodes and weights on [0,1] for the intensity integration
  !! (sum(w) = 1, flux = 2 pi sum(w mu I))
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.2113248654051871_dp, 0.7886751345948129_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)

  !! Single scattering albedo cap and minimum layer optical depth
  real(dp), parameter :: w_max = 1.0_dp - 1.0e-6_dp
  real(dp), parameter :: dtau_min = 1.0e-10_dp

  private :: lw_Fu97_method, em1x
  public :: lw_Fu97

contains

  subroutine lw_Fu97(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev, nb, ng
    real(dp), dimension(ng), intent(in) :: gw
    real(dp), dimension(nb+1), intent(in) :: wn_e
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(ng,nb,nlev), intent(in) :: tau_e
    real(dp), dimension(ng,nb,nlay), intent(in) :: ssa, gg
    real(dp), intent(in) :: Tint

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down, lw_net
    real(dp), intent(out) :: olr

    !! Work variables
    integer :: k, b, g
    real(dp), dimension(nlev) :: Te
    real(dp), dimension(nb,nlev) :: be
    real(dp), dimension(nb) :: be_int
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    !! Use WENO4 method to (smoothly) interpolate layers to levels, in ln p
    Te(:) = interpolate_weno4(log(pe), log(pl), Tl, .False.)

    !! Edges are linearly extrapolated in log(p)-log(T) space.
    Te(1) = 10.0_dp**(log10(Te(2)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * &
      & log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Te(nlay)) + (log10(pe(nlev)/pe(nlay))/ &
      & log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    !! Find integrated planck function for each band for each level
    do k = 1, nlev
      call BB_integrate(nb, Te(k), wn_e, be(:,k))
    end do
    call BB_integrate(nb, Tint, wn_e, be_int) ! Integrated planck function intensity for internal temperature

    !! Longwave flux calculation - sum over each band and g-ordinance
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp
    do b = 1, nb
      do g = 1, ng
        call lw_Fu97_method(nlay, nlev, be(b,:), be_int(b), tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), &
          & lw_up_g(:), lw_down_g(:))
        lw_up(:) = lw_up(:) + lw_up_g(:) * gw(g)
        lw_down(:) = lw_down(:) + lw_down_g(:) * gw(g)
      end do
    end do

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_Fu97

  subroutine lw_Fu97_method(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables and arrays
    integer :: k, i
    real(dp), dimension(nlay) :: dtau, w0, f, hg, g1, g2, kk, rr, aa, B1, R, T, Eu, Ed, A, C
    real(dp), dimension(nlev) :: Rt, Sd, Rb, Su, Id2, Iu2
    real(dp), dimension(nmu,nlay) :: Tq, up_src, dn_src
    real(dp), dimension(nmu,nlev) :: lw_up_g, lw_down_g
    real(dp) :: e1, e2, den, m, a1, a2, det
    real(dp) :: Ipp0, Ipm0, IppD, IpmD
    real(dp) :: Ti, pp, pm, wh, c0p, c0m, c1, alp_p, bet_p, alp_m, bet_m, Kdec, Kgro, lin1

    !! Layer optical depths and single scattering albedo
    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)
    w0(:) = min(max(w_in(:), 0.0_dp), w_max)

    !! delta-Eddington scaling, two-term phase function P = 1 + 3 g' mu mu'
    f(:) = g_in(:)**2
    dtau(:) = (1.0_dp - w0(:)*f(:))*dtau(:)
    w0(:) = w0(:)*(1.0_dp - f(:))/(1.0_dp - w0(:)*f(:))
    hg(:) = g_in(:)/(1.0_dp + g_in(:))

    !! Step 1: delta-two-stream (discrete ordinates at mu1) layer properties, intensity units:
    !! dI+/dtau = g1 I+ - g2 I- - (1-w)B/mu1,  -dI-/dtau = g1 I- - g2 I+ - (1-w)B/mu1
    do k = 1, nlay

      g1(k) = (1.0_dp - 0.5_dp*w0(k)*(1.0_dp + 3.0_dp*hg(k)*mu1**2))/mu1
      g2(k) = 0.5_dp*w0(k)*(1.0_dp - 3.0_dp*hg(k)*mu1**2)/mu1
      kk(k) = sqrt(max((g1(k) - g2(k))*(g1(k) + g2(k)), 0.0_dp))
      rr(k) = g2(k)/(g1(k) + kk(k))  ! reflection of a semi-infinite layer (mode ratio)
      aa(k) = 1.0_dp/(g1(k) + g2(k)) ! particular solution I+/- = B(tau) +/- aa B1

      if (dtau(k) <= dtau_min) then
        !! Optically negligible layer - transparent
        R(k) = 0.0_dp
        T(k) = 1.0_dp
        Eu(k) = 0.0_dp
        Ed(k) = 0.0_dp
        B1(k) = 0.0_dp
        cycle
      end if

      !! Reflection and transmission (equal from both sides)
      e1 = exp(-kk(k)*dtau(k))
      e2 = e1**2
      den = kk(k) + g1(k) + (kk(k) - g1(k))*e2
      R(k) = g2(k)*(1.0_dp - e2)/den
      T(k) = 2.0_dp*kk(k)*e1/den

      !! Thermal sources, Planck function linear in tau
      B1(k) = (be(k+1) - be(k))/dtau(k)
      Ipp0 = be(k) + aa(k)*B1(k)
      Ipm0 = be(k) - aa(k)*B1(k)
      IppD = be(k+1) + aa(k)*B1(k)
      IpmD = be(k+1) - aa(k)*B1(k)
      Eu(k) = Ipp0 - R(k)*Ipm0 - T(k)*IppD  ! emitted up at the layer top
      Ed(k) = IpmD - R(k)*IppD - T(k)*Ipm0  ! emitted down at the layer bottom

    end do

    !! Two-stream adding - composites from the top (no incident intensity)
    Rt(1) = 0.0_dp
    Sd(1) = 0.0_dp
    do k = 1, nlay
      m = 1.0_dp/(1.0_dp - Rt(k)*R(k))
      Rt(k+1) = R(k) + T(k)*m*Rt(k)*T(k)
      Sd(k+1) = Ed(k) + T(k)*m*(Sd(k) + Rt(k)*Eu(k))
    end do

    !! Composites from the bottom - gas giant: I_up = I_down + B(Tint)
    Rb(nlev) = 1.0_dp
    Su(nlev) = be_int
    do k = nlay, 1, -1
      m = 1.0_dp/(1.0_dp - Rb(k+1)*R(k))
      Rb(k) = R(k) + T(k)*m*Rb(k+1)*T(k)
      Su(k) = Eu(k) + T(k)*m*(Su(k+1) + Rb(k+1)*Ed(k))
    end do

    !! Two-stream intensities at each level
    Id2(:) = (Sd(:) + Rt(:)*Su(:))/(1.0_dp - Rt(:)*Rb(:))
    Iu2(:) = Su(:) + Rb(:)*Id2(:)

    !! Two-stream solution inside each layer (t = local depth, D = dtau):
    !! I+(t) = B(t) + aa B1 + A rr e^{-k t} + C e^{-k (D - t)}
    !! I-(t) = B(t) - aa B1 + A e^{-k t} + C rr e^{-k (D - t)}
    do k = 1, nlay
      if (dtau(k) <= dtau_min) then
        A(k) = 0.0_dp
        C(k) = 0.0_dp
        cycle
      end if
      e1 = exp(-kk(k)*dtau(k))
      a1 = Id2(k) - (be(k) - aa(k)*B1(k))
      a2 = Iu2(k+1) - (be(k+1) + aa(k)*B1(k))
      det = 1.0_dp - (rr(k)*e1)**2
      A(k) = (a1 - rr(k)*e1*a2)/det
      C(k) = (a2 - rr(k)*e1*a1)/det
    end do

    !! Step 2: four-stream source function integration along each double-Gauss angle
    do k = 1, nlay

      if (dtau(k) <= dtau_min) then
        Tq(:,k) = 1.0_dp
        up_src(:,k) = 0.0_dp
        dn_src(:,k) = 0.0_dp
        cycle
      end if

      e1 = exp(-kk(k)*dtau(k))
      wh = 0.5_dp*w0(k)

      do i = 1, nmu

        Ti = exp(-dtau(k)/uarr(i))
        Tq(i,k) = Ti

        !! Phase function between the four-stream and two-stream angles
        pp = 1.0_dp + 3.0_dp*hg(k)*uarr(i)*mu1 ! P(mu_i, mu1) = P(-mu_i, -mu1)
        pm = 1.0_dp - 3.0_dp*hg(k)*uarr(i)*mu1 ! P(mu_i,-mu1) = P(-mu_i,  mu1)

        !! S(+/-mu_i, t) = c0 + c1 t + alp e^{-k t} + bet e^{-k (D - t)}
        !! (the B(t) coefficient (1 - w) + w (pp + pm)/2 is exactly 1)
        c1 = B1(k)
        c0p = be(k) + wh*aa(k)*B1(k)*(pp - pm)
        c0m = be(k) - wh*aa(k)*B1(k)*(pp - pm)
        alp_p = wh*(pp*rr(k) + pm)*A(k)
        bet_p = wh*(pp + pm*rr(k))*C(k)
        alp_m = wh*(pm*rr(k) + pp)*A(k)
        bet_m = wh*(pm + pp*rr(k))*C(k)

        !! Kernels over the layer for the upward ray (weight e^{-t/mu}/mu); the downward ray
        !! is the mirror t -> D - t (decaying and growing exponentials swap)
        Kdec = (1.0_dp - Ti*e1)/(1.0_dp + kk(k)*uarr(i))
        Kgro = (dtau(k)/uarr(i))*exp(-min(kk(k), 1.0_dp/uarr(i))*dtau(k)) &
          & * em1x(abs(kk(k) - 1.0_dp/uarr(i))*dtau(k))
        lin1 = uarr(i)*(1.0_dp - Ti) - dtau(k)*Ti

        up_src(i,k) = c0p*(1.0_dp - Ti) + c1*lin1 + alp_p*Kdec + bet_p*Kgro
        dn_src(i,k) = (c0m + c1*dtau(k))*(1.0_dp - Ti) - c1*lin1 + alp_m*Kgro + bet_m*Kdec

      end do
    end do

    !! Four-stream sweeps: top - no incident intensity, bottom - I_up = I_down + B(Tint)
    lw_down_g(:,1) = 0.0_dp
    do k = 1, nlay
      lw_down_g(:,k+1) = lw_down_g(:,k)*Tq(:,k) + dn_src(:,k)
    end do

    lw_up_g(:,nlev) = lw_down_g(:,nlev) + be_int
    do k = nlay, 1, -1
      lw_up_g(:,k) = lw_up_g(:,k+1)*Tq(:,k) + up_src(:,k)
    end do

    !! Fluxes from the four-stream intensities
    flx_down(:) = 0.0_dp
    flx_up(:) = 0.0_dp
    do i = 1, nmu
      flx_down(:) = flx_down(:) + lw_down_g(i,:)*w(i)*uarr(i)
      flx_up(:) = flx_up(:) + lw_up_g(i,:)*w(i)*uarr(i)
    end do
    flx_down(:) = 2.0_dp*pi*flx_down(:)
    flx_up(:) = 2.0_dp*pi*flx_up(:)

  end subroutine lw_Fu97_method

  !! (1 - exp(-x))/x, stable as x -> 0
  elemental function em1x(x) result(r)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: r
    if (x < 1.0e-4_dp) then
      r = 1.0_dp - 0.5_dp*x + x*x/6.0_dp
    else
      r = (1.0_dp - exp(-x))/x
    end if
  end function em1x

end module lw_Fu97_mod