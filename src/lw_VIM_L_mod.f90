!!!
! Elspeth KH Lee - Nov 2026 : Initial version
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!
! lw: Variational Iteration Method, linear-in-tau Planck variant (alpha-4VIM_L)
!     Follows the structure of Zhang et al. (2017) but with B(t) = B_k + B1 t inside each layer.
!     Zeroth order: linear absorption approximation (AA) sweeps with
!                   eps = sqrt((1-w)(1-wg)) (delta-Eddington scaled quantities).
!     First order:  analytic scattering source from the zeroth-order intensities with a
!                   two-term phase function P = 1 + 3 g mu mu', then one final sweep with
!                   full extinction.
!
! Pros: no matrix inversions, no beta = log(B ratio) (robust in inversions and for small dtau)
! Cons: single scattering iteration on AA - degrades for w -> 1 over optically thick columns
!!!

module lw_VIM_L_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4
  use BB_integrate_mod, only : BB_integrate
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Gauss–Jacobi-5 quadrature for 4 nodes - 8 stream (Hogan 2024)
  !! (any mu-weighted quadrature with sum(w) = 1 works; weights for the scattering
  !!  integral are derived from these below)
  integer, parameter :: nmu = 4
  real(dp), dimension(nmu), parameter :: uarr = (/0.0454586727_dp, 0.2322334416_dp, 0.5740198775_dp, 0.9030775973_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.0092068785_dp, 0.1285704278_dp, 0.4323381850_dp, 0.4298845087_dp /)

  private :: lw_VIM_L_fs, em1x, q2, q3
  public :: lw_VIM_L

contains

  subroutine lw_VIM_L(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
        call lw_VIM_L_fs(nlay, nlev, be(b,:), be_int(b), tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), &
          & lw_up_g(:), lw_down_g(:))
        lw_up(:) = lw_up(:) + lw_up_g(:) * gw(g)
        lw_down(:) = lw_down(:) + lw_down_g(:) * gw(g)
      end do
    end do

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_VIM_L

  subroutine lw_VIM_L_fs(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
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
    real(dp), dimension(nlay) :: dtau, f, w0, hg, eps, dB
    real(dp), dimension(nmu) :: cw, b, y, Ap, Am, Pp, Pm, Zpp, Zpm, Hpp, Hpm
    real(dp), dimension(nmu, nlay) :: T, L, TE, G0, Sp, Sm
    real(dp), dimension(nmu, nlev) :: lw_up_g, lw_down_g
    real(dp) :: a, u, emu, q2u, q3u

    !! Calculate dtau in each layer
    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)

    !! delta-Eddington scaling - consistent with the 2-term (1 + 3 g mu mu') phase function
    f(:) = g_in(:)**2
    w0(:) = w_in(:)*(1.0_dp - f(:))/(1.0_dp - w_in(:)*f(:))
    dtau(:) = (1.0_dp - w_in(:)*f(:))*dtau(:)
    hg(:) = g_in(:)/(1.0_dp + g_in(:))

    !! Modified co-albedo for the zeroth-order AA solution (used consistently everywhere below)
    eps(:) = sqrt((1.0_dp - w0(:))*(1.0_dp - w0(:)*hg(:)))

    !! Planck difference across each layer (B1 * dtau)
    dB(:) = be(2:nlev) - be(1:nlay)

    !! Half-range d(mu) weights for the scattering integral from the mu-weighted flux weights
    cw(:) = w(:)/uarr(:)
    cw(:) = cw(:)/sum(cw(:))

    !! Transmission and linear-source terms
    do k = 1, nlay
      do i = 1, nmu
        TE(i,k) = exp(-eps(k)*dtau(k)/uarr(i))       ! AA transmission
        G0(i,k) = dB(k)*em1x(eps(k)*dtau(k)/uarr(i)) ! = mu B1 (1 - TE)/eps
        T(i,k) = exp(-dtau(k)/uarr(i))               ! full extinction transmission
        L(i,k) = dB(k)*em1x(dtau(k)/uarr(i))         ! = mu B1 (1 - T)
      end do
    end do

    !! Zeroth order - linear AA sweeps
    ! Top boundary condition - 0 flux downward from top boundary
    lw_down_g(:,1) = 0.0_dp
    do k = 1, nlay
      lw_down_g(:,k+1) = lw_down_g(:,k)*TE(:,k) + (be(k+1) - be(k)*TE(:,k)) - G0(:,k)
    end do

    ! Lower boundary condition - internal heat definition Fint = F_up - F_down (intensity units)
    lw_up_g(:,nlev) = lw_down_g(:,nlev) + be_int
    do k = nlay, 1, -1
      lw_up_g(:,k) = lw_up_g(:,k+1)*TE(:,k) + (be(k) - be(k+1)*TE(:,k)) + G0(:,k)
    end do

    !! If no scattering component in profile, AA (eps = 1) is the exact linear solution - return
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
    ! Zeroth-order intensity inside layer (t = local depth, D = dtau):
    !  I0+_j(t) = Ap_j e^{-eps(D-t)/mu_j} + B(t) + (mu_j B1/eps)(1 - e^{-eps(D-t)/mu_j})
    !  I0-_j(t) = Am_j e^{-eps t/mu_j}    + B(t) - (mu_j B1/eps)(1 - e^{-eps t/mu_j})
    ! Kernels (upward i, weight e^{-t/mu_i}/mu_i; downward i by symmetry t -> D - t):
    !  Zpp = int e^{-eps(D-t)/mu_j} e^{-t/mu_i} dt/mu_i
    !  Zpm = int e^{-eps t/mu_j}    e^{-t/mu_i} dt/mu_i
    !  Hpp = (mu_j/D) int e^{-t/mu_i}(1 - e^{-eps(D-t)/mu_j})/eps dt/mu_i = ((1-T_i) - Zpp)/y
    !  Hpm = (mu_j/D) int e^{-t/mu_i}(1 - e^{-eps t/mu_j})/eps    dt/mu_i = ((1-T_i) - Zpm)/y
    ! with y = eps D/mu_j (series in y for y -> 0)
    Sp(:,:) = 0.0_dp
    Sm(:,:) = 0.0_dp

    do k = 1, nlay

      ! Cycle if no scattering component
      if (w0(k) <= 1.0e-6_dp) then
        cycle
      end if

      Ap(:) = lw_up_g(:,k+1) - be(k+1)  ! homogeneous amplitude at layer base
      Am(:) = lw_down_g(:,k) - be(k)    ! homogeneous amplitude at layer top
      b(:) = eps(k)/uarr(:)
      y(:) = eps(k)*dtau(k)/uarr(:)

      do i = 1, nmu

        a = 1.0_dp/uarr(i)
        u = dtau(k)*a
        emu = em1x(u)
        q2u = q2(u)
        q3u = q3(u)

        do j = 1, nmu
          Zpp(j) = u * exp(-min(a,b(j))*dtau(k)) * em1x(abs(a - b(j))*dtau(k))
          Zpm(j) = u * em1x((a + b(j))*dtau(k))
          if (y(j) < 1.0e-4_dp) then
            Hpp(j) = (1.0_dp - emu) - y(j)*q2u
            Hpm(j) = (emu - T(i,k)) - y(j)*q3u
          else
            Hpp(j) = ((1.0_dp - T(i,k)) - Zpp(j))/y(j)
            Hpm(j) = ((1.0_dp - T(i,k)) - Zpm(j))/y(j)
          end if
        end do

        Pp(:) = 1.0_dp + 3.0_dp*hg(k)*uarr(i)*uarr(:) ! same hemisphere
        Pm(:) = 1.0_dp - 3.0_dp*hg(k)*uarr(i)*uarr(:) ! opposite hemisphere

        ! Scattering source minus the w*Pi part (absorbed into the final sweep Planck term)
        Sp(i,k) = 0.5_dp*w0(k) * sum(cw(:) * ( Pp(:)*(Ap(:)*Zpp(:) + dB(k)*Hpp(:)) &
          & + Pm(:)*(Am(:)*Zpm(:) - dB(k)*Hpm(:)) ))
        Sm(i,k) = 0.5_dp*w0(k) * sum(cw(:) * ( Pp(:)*(Am(:)*Zpp(:) - dB(k)*Hpp(:)) &
          & + Pm(:)*(Ap(:)*Zpm(:) + dB(k)*Hpm(:)) ))

      end do
    end do

    !! Final sweeps - full extinction, (1-w)Pi + w Pi = Pi plus scattering source
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp

    lw_down_g(:,1) = 0.0_dp
    do k = 1, nlay
      lw_down_g(:,k+1) = lw_down_g(:,k)*T(:,k) + (be(k+1) - be(k)*T(:,k)) - L(:,k) + Sm(:,k)
    end do

    lw_up_g(:,nlev) = lw_down_g(:,nlev) + be_int
    do k = nlay, 1, -1
      lw_up_g(:,k) = lw_up_g(:,k+1)*T(:,k) + (be(k) - be(k+1)*T(:,k)) + L(:,k) + Sp(:,k)
    end do

    do i = 1, nmu
      flx_down(:) = flx_down(:) + lw_down_g(i,:) * w(i)
      flx_up(:) = flx_up(:) + lw_up_g(i,:) * w(i)
    end do

    !! The flux is the integrated intensity * pi (in this GJ weighting scheme)
    flx_down(:) = pi * flx_down(:)
    flx_up(:) = pi * flx_up(:)

  end subroutine lw_VIM_L_fs

  !! (1 - exp(-x))/x, stable as x -> 0
  elemental function em1x(x) result(r)
    real(dp), intent(in) :: x
    real(dp) :: r
    if (x < 1.0e-4_dp) then
      r = 1.0_dp - 0.5_dp*x + x*x/6.0_dp
    else
      r = (1.0_dp - exp(-x))/x
    end if
  end function em1x

  !! 1/2 - (1 - em1x(u))/u  (O(y) term of Hpp)
  elemental function q2(u) result(r)
    real(dp), intent(in) :: u
    real(dp) :: r
    if (u < 1.0e-3_dp) then
      r = u/6.0_dp - u*u/24.0_dp + u**3/120.0_dp
    else
      r = 0.5_dp - (1.0_dp - em1x(u))/u
    end if
  end function q2

  !! (em1x(u) - e^{-u})/u - e^{-u}/2  (O(y) term of Hpm)
  elemental function q3(u) result(r)
    real(dp), intent(in) :: u
    real(dp) :: r
    if (u < 1.0e-3_dp) then
      r = u/6.0_dp - u*u/8.0_dp + u**3/20.0_dp
    else
      r = (em1x(u) - exp(-u))/u - 0.5_dp*exp(-u)
    end if
  end function q3

end module lw_VIM_L_mod