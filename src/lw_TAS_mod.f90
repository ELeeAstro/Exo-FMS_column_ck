!!!
! Elspeth KH Lee - Nov 2026 : Initial version
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!                - Nov 2026 : Optimisations - backscatter series coefficients once per call, series
!                             stopped once g^l is negligible, (1 - T)/x source term from the
!                             transmission, contiguous per g-ordinance inputs
!
! lw: Chou et al. (1999) scaling with the Tang et al. (2018) adjustment scheme (TAS),
!     thermal emission only. Used for longwave scattering in RTTOV and as the LW rescaling
!     option in RTE+RRTMGP.
!     Chou scaling: the diffuse field is replaced by I(mu') = I(mu) in the same hemisphere and
!       B in the opposite hemisphere, so the scattering only removes the backscattered fraction b:
!         mu/(1 - w(1-b)) dI/dtau = I - B,  i.e. a no-scattering problem in the scaled optical
!         depth tau' = (1 - w(1-b)) tau = tau_abs + b tau_sca.
!     TAS: for the upward intensity the opposite-hemisphere B is replaced by the (Chou) downward
!       intensity at the same |mu|:
!         mu dI+/dtau' = I+ - B - kappa (I-(tau',mu) - B),  kappa = w b/(1 - w(1-b)).
!       With B linear in tau' and the Chou downward intensity inside the layer, the correction
!       integrates analytically to (t = 0 layer top, D' = scaled layer optical depth, T = e^{-D'/mu}):
!         dI+(0) = f_T kappa [ (I-(0) - B_top) - (I-(D') - B_bot) T - mu B1' (1 - T) ],
!       where f_T = 1/2 is the exact solution of the TAS equations, which is the default here.
!       Tang et al. (2018) suggested a tunable "Tang factor" (0.3), better for weak scattering only.
!     Downward intensities are the Chou solution (no correction), as in Tang et al. (2018).
!
! Notes:
!   - b is the hemispheric backscatter fraction of the full (unscaled) Henyey-Greenstein phase
!     function, b = 1/2 int_0^1 dmu int_-1^0 Pbar(mu,mu') dmu' = 1/2 - 1/2 sum_odd (2l+1) g^l c_l^2,
!     with c_l = int_0^1 P_l dmu = (P_(l-1)(0) - P_(l+1)(0))/(2l+1); no delta-scaling is applied.
!   - tang_fact = 0 gives the plain Chou et al. (1999) scaling.
!   - Hogan (2024) Gauss-Jacobi quadrature for the angular integration (flux = pi sum(w I)).
!
! Pros: no-scattering cost (two sweeps, no matrix solves), much better than Chou scaling alone
! Cons: single-pass correction - fails for strongly scattering optically thick atmospheres (w -> 1);
!       the Tang factor is empirical
!!!

module lw_TAS_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4
  use BB_integrate_mod, only : BB_integrate
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Tang factor: 0.5 = exact TAS solution, 0.3 = Tang et al. (2018) suggestion, 0 = Chou scaling only
  real(dp), parameter :: tang_fact = 0.5_dp

  !! Number of Legendre terms for the HG backscatter fraction series, and |g|^l below which the
  !! series is stopped (the remaining terms are < 1e-17 for |g| < 0.86; larger |g| use all nleg terms)
  integer, parameter :: nleg = 256
  real(dp), parameter :: gl_tol = 1.0e-17_dp

  !! Gauss–Jacobi-5 quadrature for 4 nodes - 8 stream (Hogan 2024)
  integer, parameter :: nmu = 4
  real(dp), dimension(nmu), parameter :: uarr = (/0.0454586727_dp, 0.2322334416_dp, 0.5740198775_dp, 0.9030775973_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.0092068785_dp, 0.1285704278_dp, 0.4323381850_dp, 0.4298845087_dp /)
  real(dp), dimension(nmu), parameter :: ruarr = 1.0_dp/uarr

  private :: lw_TAS_method, em1x
  public :: lw_TAS

contains

  subroutine lw_TAS(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
    integer :: k, b, g, l
    real(dp), dimension(nlev) :: Te
    real(dp), dimension(nb,nlev) :: be
    real(dp), dimension(nb) :: be_int
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g
    real(dp), dimension(nlev) :: be_b, tau_g
    real(dp), dimension(nlay) :: ssa_g, gg_g
    real(dp), dimension(0:nleg+1) :: P0
    real(dp), dimension(nleg) :: cl2

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

    !! (2l+1) c_l^2 for odd l, c_l = int_0^1 P_l dmu = (P_(l-1)(0) - P_(l+1)(0))/(2l+1)
    P0(0) = 1.0_dp
    P0(1) = 0.0_dp
    do l = 1, nleg
      P0(l+1) = -real(l,dp)/real(l+1,dp)*P0(l-1)
    end do
    cl2(:) = 0.0_dp
    do l = 1, nleg, 2
      cl2(l) = (P0(l-1) - P0(l+1))**2/real(2*l+1,dp)
    end do

    !! Longwave flux calculation - sum over each band and g-ordinance
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp
    do b = 1, nb
      be_b(:) = be(b,:)
      do g = 1, ng
        ! Contiguous copies (avoids array temporaries for the strided slices)
        tau_g(:) = tau_e(g,b,:)
        ssa_g(:) = ssa(g,b,:)
        gg_g(:) = gg(g,b,:)
        call lw_TAS_method(nlay, nlev, be_b, be_int(b), tau_g, ssa_g, gg_g, cl2, lw_up_g(:), lw_down_g(:))
        lw_up(:) = lw_up(:) + lw_up_g(:) * gw(g)
        lw_down(:) = lw_down(:) + lw_down_g(:) * gw(g)
      end do
    end do

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_TAS

  subroutine lw_TAS_method(nlay, nlev, be, be_int, tau_in, w_in, g_in, cl2, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int
    real(dp), dimension(nleg), intent(in) :: cl2   ! (2l+1) c_l^2 (odd l)

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables and arrays
    integer :: k, i, l
    real(dp), dimension(nlay) :: w0, b, sc, dtau, kap, dB
    real(dp), dimension(nmu,nlay) :: T, Ls
    real(dp), dimension(nmu,nlev) :: lw_up_g, lw_down_g
    real(dp), dimension(nmu) :: dm0, dmD
    real(dp) :: s, gl, g2, x

    !! Single scattering albedo and HG hemispheric backscatter fraction
    w0(:) = min(max(w_in(:), 0.0_dp), 1.0_dp)
    do k = 1, nlay
      s = 0.0_dp
      gl = g_in(k)
      g2 = g_in(k)**2
      do l = 1, nleg, 2
        if (abs(gl) < gl_tol) exit
        s = s + cl2(l)*gl
        gl = gl*g2
      end do
      b(k) = min(max(0.5_dp - 0.5_dp*s, 0.0_dp), 1.0_dp)
    end do

    !! Chou scaling factor, scaled optical depth and the TAS coefficient kappa
    sc(:) = max(1.0_dp - w0(:)*(1.0_dp - b(:)), 1.0e-12_dp)
    dtau(:) = sc(:)*(tau_in(2:nlev) - tau_in(1:nlay))
    kap(:) = w0(:)*b(:)/sc(:)

    !! Transmission and linear-in-tau' source terms
    dB(:) = be(2:nlev) - be(1:nlay)
    do k = 1, nlay
      do i = 1, nmu
        x = dtau(k)*ruarr(i)
        T(i,k) = exp(-x)
        ! = mu B1' (1 - T) = dB (1 - T)/x
        if (x < 1.0e-4_dp) then
          Ls(i,k) = dB(k)*em1x(x)
        else
          Ls(i,k) = dB(k)*(1.0_dp - T(i,k))/x
        end if
      end do
    end do

    !! Downward sweep - Chou scaling (no-scattering solution in tau')
    ! Top boundary condition - no incident intensity
    lw_down_g(:,1) = 0.0_dp
    do k = 1, nlay
      lw_down_g(:,k+1) = lw_down_g(:,k)*T(:,k) + (be(k+1) - be(k)*T(:,k)) - Ls(:,k)
    end do

    !! Upward sweep - Chou scaling plus the Tang adjustment term
    ! Lower boundary condition - internal heat definition Fint = F_up - F_down
    lw_up_g(:,nlev) = lw_down_g(:,nlev) + be_int
    do k = nlay, 1, -1
      dm0(:) = lw_down_g(:,k) - be(k)     ! (I- - B) at the layer top
      dmD(:) = lw_down_g(:,k+1) - be(k+1) ! (I- - B) at the layer bottom
      lw_up_g(:,k) = lw_up_g(:,k+1)*T(:,k) + (be(k) - be(k+1)*T(:,k)) + Ls(:,k) &
        & + tang_fact*kap(k)*(dm0(:) - dmD(:)*T(:,k) - Ls(:,k))
    end do

    !! Fluxes (Gauss-Jacobi weighting: flux = pi * sum(w I))
    flx_down(:) = 0.0_dp
    flx_up(:) = 0.0_dp
    do i = 1, nmu
      flx_down(:) = flx_down(:) + lw_down_g(i,:)*w(i)
      flx_up(:) = flx_up(:) + lw_up_g(i,:)*w(i)
    end do
    flx_down(:) = pi*flx_down(:)
    flx_up(:) = pi*flx_up(:)

  end subroutine lw_TAS_method

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

end module lw_TAS_mod