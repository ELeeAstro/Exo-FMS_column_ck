!!!
! Elspeth KH Lee - Aug 2023 : Initial version
!                - Nov 2026 : Renamed lw_VIM_E; phase function pairing and scattering-integral
!                             weights corrected; rewritten in the style of lw_VIM_L
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!
! lw: Variational Iteration Method, exponential-in-tau Planck variant (alpha-4VIM_E)
!     Follows Zhang et al. (2017, J. Atmos. Sci. 74, 419) with B(t) = B_k e^{-beta t} inside
!     each layer, beta = -ln(B_{k+1}/B_k)/dtau.
!     Zeroth order: AA sweeps with the extended AA co-albedo eps_g = sqrt((1-w)(1-wg))
!                   (Li 2002) for the level intensities.
!     First order:  analytic scattering source from the zeroth-order intensities with a
!                   two-term phase function P = 1 + 3 g mu mu', then one final sweep with
!                   full extinction.
!
! Pros: no matrix inversions, OLR much better than AA alone
! Cons: single scattering iteration on AA - degrades for w -> 1 over optically thick columns
!!!

module lw_VIM_E_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4
  use BB_integrate_mod, only : BB_integrate
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Minimum layer optical depth (beta = -ln(B ratio)/dtau)
  real(dp), parameter :: dtau_floor = 1.0e-12_dp

  !! c+/- resonance guard: |eps -/+ mu beta| is held no closer to zero than res_floor*eps.
  !! Results are only weakly sensitive to the value (0.1 and 0.3 differ by ~3% in the
  !! worst-case heating-rate error).
  real(dp), parameter :: res_floor = 0.1_dp

  !! Gauss–Jacobi-5 quadrature for 4 nodes - 8 stream (Hogan 2024)
  !! (any mu-weighted quadrature with sum(w) = 1 works; weights for the scattering
  !!  integral are derived from these below)
  integer, parameter :: nmu = 4
  real(dp), dimension(nmu), parameter :: uarr = (/0.0454586727_dp, 0.2322334416_dp, 0.5740198775_dp, 0.9030775973_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.0092068785_dp, 0.1285704278_dp, 0.4323381850_dp, 0.4298845087_dp /)

  private :: lw_VIM_E_fs, aa_source_down, aa_source_up, em1x, exprel
  public :: lw_VIM_E

contains

  subroutine lw_VIM_E(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
        call lw_VIM_E_fs(nlay, nlev, be(b,:), be_int(b), tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), &
          & lw_up_g(:), lw_down_g(:))
        lw_up(:) = lw_up(:) + lw_up_g(:) * gw(g)
        lw_down(:) = lw_down(:) + lw_down_g(:) * gw(g)
      end do
    end do

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_VIM_E

  subroutine lw_VIM_E_fs(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables and arrays
    integer :: k, i, j
    real(dp), dimension(nlay) :: dtau, f, w0, hg, epsg, eps, beta
    real(dp), dimension(nmu) :: cw, b, cp, cm, Ap, Am, Pp, Pm, Zpp, Zpm
    real(dp), dimension(nmu, nlay) :: TG, TE, T, Sp, Sm
    real(dp), dimension(nmu, nlev) :: lw_up_g, lw_down_g
    real(dp) :: a, u, d, Pi_up, Pi_dn

    !! Calculate dtau in each layer
    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)

    !! delta-Eddington scaling - consistent with the 2-term (1 + 3 g mu mu') phase function
    f(:) = g_in(:)**2
    w0(:) = w_in(:)*(1.0_dp - f(:))/(1.0_dp - w_in(:)*f(:))
    dtau(:) = (1.0_dp - w_in(:)*f(:))*dtau(:)
    hg(:) = g_in(:)/(1.0_dp + g_in(:))

    !! Avoid singular beta in numerically transparent layers (after scaling: the transport depth)
    dtau(:) = max(dtau(:), dtau_floor)

    !! Exponential-in-tau Planck function: B(t) = B_k e^{-beta t}
    beta(:) = -log(be(2:nlev)/be(1:nlay))/dtau(:)

    !! Co-albedos: extended AA for the zeroth-order level intensities, plain AA in-layer profile
    epsg(:) = sqrt((1.0_dp - w0(:))*(1.0_dp - w0(:)*hg(:)))
    eps(:) = 1.0_dp - w0(:)

    !! Half-range d(mu) weights for the scattering integral from the mu-weighted flux weights
    cw(:) = w(:)/uarr(:)
    cw(:) = cw(:)/sum(cw(:))

    !! Transmission functions
    do k = 1, nlay
      do i = 1, nmu
        TG(i,k) = exp(-epsg(k)*dtau(k)/uarr(i)) ! zeroth-order (extended AA) transmission
        TE(i,k) = exp(-eps(k)*dtau(k)/uarr(i))  ! in-layer AA transmission
        T(i,k) = exp(-dtau(k)/uarr(i))          ! full extinction transmission
      end do
    end do

    !! Zeroth order - extended AA sweeps
    ! Top boundary condition - 0 flux downward from top boundary
    lw_down_g(:,1) = 0.0_dp
    do k = 1, nlay
      lw_down_g(:,k+1) = lw_down_g(:,k)*TG(:,k) + &
        & aa_source_down(be(k), be(k+1), epsg(k)*dtau(k)/uarr(:), TG(:,k))
    end do

    ! Lower boundary condition - internal heat definition Fint = F_up - F_down (intensity units)
    lw_up_g(:,nlev) = lw_down_g(:,nlev) + be_int
    do k = nlay, 1, -1
      lw_up_g(:,k) = lw_up_g(:,k+1)*TG(:,k) + &
        & aa_source_up(be(k), be(k+1), epsg(k)*dtau(k)/uarr(:), TG(:,k))
    end do

    !! If no scattering component in profile, the zeroth order is the solution - return
    if (all(w0(:) <= 1.0e-6_dp)) then
      flx_down(:) = 0.0_dp
      flx_up(:) = 0.0_dp
      do i = 1, nmu
        flx_down(:) = flx_down(:) + lw_down_g(i,:) * w(i)
        flx_up(:) = flx_up(:) + lw_up_g(i,:) * w(i)
      end do
      flx_down(:) = pi * flx_down(:)
      flx_up(:) = pi * flx_up(:)
      return
    end if

    !! First order - VIM scattering source functions
    ! Zeroth-order intensity inside layer (t = local depth, D = dtau, AA with eps = 1 - w):
    !  I0+_j(t) = Ap_j e^{-eps(D-t)/mu_j} + c+_j B(t),  c+_j = eps/(eps + mu_j beta)
    !  I0-_j(t) = Am_j e^{-eps t/mu_j}    + c-_j B(t),  c-_j = eps/(eps - mu_j beta)
    ! with Ap_j = I0+_j(D) - c+_j B_{k+1} and Am_j = I0-_j(0) - c-_j B_k from the level sweeps.
    ! Kernels (upward i, weight e^{-t/mu_i}/mu_i; downward i by symmetry t -> D - t):
    !  Zpp   = int e^{-eps(D-t)/mu_j} e^{-t/mu_i} dt/mu_i
    !  Zpm   = int e^{-eps t/mu_j}    e^{-t/mu_i} dt/mu_i
    !  Pi_up = int B(t) e^{-t/mu_i} dt/mu_i,  Pi_dn = int B(t) e^{-(D-t)/mu_i} dt/mu_i
    Sp(:,:) = 0.0_dp
    Sm(:,:) = 0.0_dp

    do k = 1, nlay

      ! Cycle if no scattering component
      if (w0(k) <= 1.0e-6_dp) then
        cycle
      end if

      ! Particular-solution amplitudes, held off their resonances
      do j = 1, nmu
        d = eps(k) + uarr(j)*beta(k)
        if (abs(d) < res_floor*eps(k)) d = sign(res_floor*eps(k), d)
        cp(j) = eps(k)/d
        d = eps(k) - uarr(j)*beta(k)
        if (abs(d) < res_floor*eps(k)) d = sign(res_floor*eps(k), d)
        cm(j) = eps(k)/d
      end do

      Ap(:) = lw_up_g(:,k+1) - cp(:)*be(k+1)  ! homogeneous amplitude at layer base
      Am(:) = lw_down_g(:,k) - cm(:)*be(k)    ! homogeneous amplitude at layer top
      b(:) = eps(k)/uarr(:)

      do i = 1, nmu

        a = 1.0_dp/uarr(i)
        u = dtau(k)*a

        do j = 1, nmu
          Zpp(j) = u * exp(-min(a,b(j))*dtau(k)) * em1x(abs(a - b(j))*dtau(k))
          Zpm(j) = u * em1x((a + b(j))*dtau(k))
        end do

        Pi_up = be(k)*u*em1x((a + beta(k))*dtau(k))
        Pi_dn = be(k+1)*u*em1x((a - beta(k))*dtau(k))

        Pp(:) = 1.0_dp + 3.0_dp*hg(k)*uarr(i)*uarr(:) ! same hemisphere
        Pm(:) = 1.0_dp - 3.0_dp*hg(k)*uarr(i)*uarr(:) ! opposite hemisphere

        Sp(i,k) = 0.5_dp*w0(k) * sum(cw(:) * ( Pp(:)*(Ap(:)*Zpp(:) + cp(:)*Pi_up) &
          & + Pm(:)*(Am(:)*Zpm(:) + cm(:)*Pi_up) ))
        Sm(i,k) = 0.5_dp*w0(k) * sum(cw(:) * ( Pp(:)*(Am(:)*Zpp(:) + cm(:)*Pi_dn) &
          & + Pm(:)*(Ap(:)*Zpm(:) + cp(:)*Pi_dn) ))

      end do
    end do

    !! Final sweeps - full extinction, (1-w) Planck emission plus scattering source
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp

    lw_down_g(:,1) = 0.0_dp
    do k = 1, nlay
      lw_down_g(:,k+1) = lw_down_g(:,k)*T(:,k) + &
        & eps(k)*aa_source_down(be(k), be(k+1), dtau(k)/uarr(:), T(:,k)) + Sm(:,k)
    end do

    lw_up_g(:,nlev) = lw_down_g(:,nlev) + be_int
    do k = nlay, 1, -1
      lw_up_g(:,k) = lw_up_g(:,k+1)*T(:,k) + &
        & eps(k)*aa_source_up(be(k), be(k+1), dtau(k)/uarr(:), T(:,k)) + Sp(:,k)
    end do

    do i = 1, nmu
      flx_down(:) = flx_down(:) + lw_down_g(i,:) * w(i)
      flx_up(:) = flx_up(:) + lw_up_g(i,:) * w(i)
    end do

    !! The flux is the integrated intensity * pi (in this GJ weighting scheme)
    flx_down(:) = pi * flx_down(:)
    flx_up(:) = pi * flx_up(:)

  end subroutine lw_VIM_E_fs

  !! int_0^D B(t) e^{-x(D-t)/D} x dt/D for B exponential in t between btop and bbot
  !! (downward emission of a layer of slant optical depth x, trans = e^{-x})
  elemental function aa_source_down(btop, bbot, x, trans) result(source)
    real(dp), intent(in) :: btop, bbot, x, trans
    real(dp) :: source
    real(dp) :: z
    z = x + log(bbot/btop)
    if (abs(z) < 1.0e-4_dp) then
      source = btop*trans*x*exprel(z)
    else
      source = (x/z)*(bbot - btop*trans)
    end if
  end function aa_source_down

  !! As aa_source_down for the upward emission
  elemental function aa_source_up(btop, bbot, x, trans) result(source)
    real(dp), intent(in) :: btop, bbot, x, trans
    real(dp) :: source
    real(dp) :: z
    z = log(bbot/btop) - x
    if (abs(z) < 1.0e-4_dp) then
      source = btop*x*exprel(z)
    else
      source = (x/(x - log(bbot/btop)))*(btop - bbot*trans)
    end if
  end function aa_source_up

  !! (1 - exp(-x))/x, stable as x -> 0 (either sign)
  elemental function em1x(x) result(r)
    real(dp), intent(in) :: x
    real(dp) :: r
    if (abs(x) < 1.0e-4_dp) then
      r = 1.0_dp - x*(0.5_dp - x*(1.0_dp/6.0_dp - x/24.0_dp))
    else
      r = (1.0_dp - exp(-x))/x
    end if
  end function em1x

  !! (exp(x) - 1)/x, stable as x -> 0
  elemental function exprel(x) result(r)
    real(dp), intent(in) :: x
    real(dp) :: r
    if (abs(x) < 1.0e-4_dp) then
      r = 1.0_dp + x*(0.5_dp + x*(1.0_dp/6.0_dp + x*(1.0_dp/24.0_dp + x/120.0_dp)))
    else
      r = (exp(x) - 1.0_dp)/x
    end if
  end function exprel

end module lw_VIM_E_mod
