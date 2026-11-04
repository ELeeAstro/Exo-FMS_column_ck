!!!
! Elspeth KH Lee - Jun 2024 : Initial version
!                - Nov 2026 : Rewritten following PICASO exactly, enabled
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!                - Nov 2026 : Optimisations - compile-time band widths in the solver, contiguous per
!                             g-ordinance inputs
!
! sw: Four-stream spherical harmonic method (SH4) with multiple scattering (Rooney et al. 2024)
!     Taken from PICASO (Batalha et al. 2019; Rooney et al. 2024) and converted to Fortran:
!     fluxes.get_reflected_SH and fluxes.setup_4_stream_fluxes, checked against PICASO commit 4524b70
!     (Oct 2026) to 1e-14. Delta-M scaling (f = g**4) and one-term HG phase function as PICASO.
!     Differences to PICASO: no Rayleigh, gas giant lower boundary, own banded solver,
!     fluxes found row by row, single zenith angle (no zcorr). None change the solution.
!     TTHG option: as PICASO, but delta-M uses the TTHG chi_nstr for w0 and dtau too (consistent).
!     Pros: Most accurate sw scheme here, works well for clouds and hazes
!     Cons: Banded matrix solve, ~5-8x slower than the two-stream methods
!!!

module sw_SH4_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Single scattering albedo cap - exactly conservative scattering (w = 1, eigenvalue k = 0) is singular.
  !! 1 - 1e-12 absorbs < 1e-8 of the incident flux even in a tau ~ 6000 column
  real(dp), parameter :: w_max = 1.0_dp - 1.0e-12_dp

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: fourpi = 4.0_dp * pi

  !! Number of streams
  integer, parameter :: nstr = 4

  !! Use Two-Term HG function (TTHG), parameters are the PICASO defaults
  !! forward fraction = frac_a + frac_b*g_back**frac_c, g_forward = c_forward*g, g_back = c_back*g
  logical, parameter :: TTHG = .False.
  real(dp), parameter :: tthg_forward = 1.0_dp, tthg_back = -0.5_dp
  real(dp), parameter :: tthg_frac_a = 1.0_dp, tthg_frac_b = -1.0_dp
  integer, parameter :: tthg_frac_c = 2

  !! Band widths of the four-stream matrix
  integer, parameter :: kl = 5, ku = 5

  private :: sw_SH_four_stream, band_solve
  public :: sw_SH4

contains

  subroutine sw_SH4(nlay, nlev, nb, ng, gw, tau_e, mu_z, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
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
          call sw_SH_four_stream(nlay, nlev, Finc(b), mu_z(nlev), tau_g, ssa_g, gg_g, &
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

  end subroutine sw_SH4

  subroutine sw_SH_four_stream(nlay, nlev, F0PI, u0, tau_in, w_in, g_in, flx_down, flx_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0PI, u0
    real(dp), dimension(nlev), intent(in) :: tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    !! Work variables
    integer :: k, l, m, n
    real(dp) :: Del, b_top, b_surface, b_surface_SH4
    real(dp), dimension(nstr) :: Pu0
    real(dp), dimension(nlay) :: dtau, w0, f_deltaM, g_f, g_b, f_f
    real(dp), dimension(nstr,nlay) :: chi
    real(dp), dimension(nstr,nlay) :: w_multi, a, bb, Dels, eta
    real(dp), dimension(nlay) :: beta, gama, disc, lam1, lam2, exptrm1, exptrm2
    real(dp), dimension(nlay) :: R1, R2, Q1, Q2, S1, S2
    real(dp), dimension(nlay) :: p1pl, p2pl, q1pl, q2pl, p1mn, p2mn, q1mn, q2mn
    real(dp), dimension(nlay) :: z1pl, z1mn, z2pl, z2mn
    real(dp), dimension(nlay) :: f00, f01, f02, f03, f10, f11, f12, f13
    real(dp), dimension(nlay) :: f20, f21, f22, f23, f30, f31, f32, f33
    real(dp), dimension(nlay) :: z1mn_up, z2mn_up, z1pl_up, z2pl_up
    real(dp), dimension(nlay) :: z1mn_down, z2mn_down, z1pl_down, z2pl_down
    real(dp), dimension(nlev) :: tau, expon, cum_trans
    real(dp), dimension(2*kl+ku+1, 4*nlay) :: Mb
    real(dp), dimension(4*nlay) :: B, X

    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w_in(:) <= 1.0e-12_dp)) then
      cum_trans(:) = tau_in(:)/u0
      flx_down(:) = F0PI * u0 * exp(-cum_trans(:))
      flx_up(:) = 0.0_dp
      return
    end if

    n = 4*nlay

    !! Phase function Legendre moments chi_l (OTHG: chi_l = g**l, TTHG as PICASO)
    if (TTHG .eqv. .False.) then
      do l = 1, nstr
        chi(l,:) = g_in(:)**l
      end do
    else
      g_f(:) = tthg_forward*g_in(:)
      g_b(:) = tthg_back*g_in(:)
      f_f(:) = tthg_frac_a + tthg_frac_b*g_b(:)**tthg_frac_c
      do l = 1, nstr
        chi(l,:) = f_f(:)*g_f(:)**l + (1.0_dp - f_f(:))*g_b(:)**l
      end do
    end if

    !! Delta-M with f = chi_nstr (g**nstr for OTHG as PICASO, the TTHG moment for TTHG)
    f_deltaM(:) = chi(nstr,:)
    w0(:) = w_in(:)*(1.0_dp - f_deltaM(:))/(1.0_dp - w_in(:)*f_deltaM(:))
    w0(:) = min(w0(:), w_max)
    dtau(:) = (tau_in(2:nlev) - tau_in(1:nlay))*(1.0_dp - w_in(:)*f_deltaM(:))
    tau(1) = 0.0_dp
    do k = 1, nlay
      tau(k+1) = tau(k) + dtau(k)
    end do

    !! Legendre polynomials at -u0
    Pu0(1) = 1.0_dp
    Pu0(2) = -u0
    Pu0(3) = (3.0_dp*u0**2 - 1.0_dp)/2.0_dp
    Pu0(4) = (5.0_dp*(-u0)**3 - 3.0_dp*(-u0))/2.0_dp

    !! Phase function moments with delta-M
    w_multi(1,:) = 1.0_dp
    do l = 1, nstr-1
      w_multi(l+1,:) = (real(2*l+1,dp)*chi(l,:) - real(2*l+1,dp)*f_deltaM(:))/(1.0_dp - f_deltaM(:))
    end do

    do l = 0, nstr-1
      a(l+1,:) = real(2*l+1,dp) - w0(:)*w_multi(l+1,:)
      bb(l+1,:) = (F0PI*(w0(:)*w_multi(l+1,:)))*Pu0(l+1)/fourpi
    end do

    !! Boundary conditions - nothing enters at the top, gas giant lower boundary
    b_top = 0.0_dp
    b_surface = 0.0_dp
    b_surface_SH4 = 0.0_dp

    !! Eigenvalues (lam2 from lam1**2 lam2**2 = gama to avoid cancellation)
    beta(:) = a(1,:)*a(2,:) + 4.0_dp*a(1,:)*a(4,:)/9.0_dp + a(3,:)*a(4,:)/9.0_dp
    gama(:) = a(1,:)*a(2,:)*a(3,:)*a(4,:)/9.0_dp
    disc(:) = sqrt(max(beta(:)**2 - 4.0_dp*gama(:), 0.0_dp))
    lam1(:) = sqrt((beta(:) + disc(:))/2.0_dp)
    lam2(:) = sqrt(2.0_dp*gama(:)/(beta(:) + disc(:)))

    !! Direct beam particular solution
    do k = 1, nlay
      Del = 9.0_dp*((1.0_dp/u0)**4 - beta(k)*(1.0_dp/u0)**2 + gama(k))
      Dels(1,k) = (a(2,k)*bb(1,k) - bb(2,k)/u0)*(a(3,k)*a(4,k) - 9.0_dp/u0**2) &
        & + 2.0_dp*(a(4,k)*bb(3,k) - 2.0_dp*a(4,k)*bb(1,k) - 3.0_dp*bb(4,k)/u0)/u0**2
      Dels(2,k) = (a(1,k)*bb(2,k) - bb(1,k)/u0)*(a(3,k)*a(4,k) - 9.0_dp/u0**2) &
        & - 2.0_dp*a(1,k)*(a(4,k)*bb(3,k) - 3.0_dp*bb(4,k)/u0)/u0
      Dels(3,k) = (a(4,k)*bb(3,k) - 3.0_dp*bb(4,k)/u0)*(a(1,k)*a(2,k) - 1.0_dp/u0**2) &
        & - 2.0_dp*a(4,k)*(a(1,k)*bb(2,k) - bb(1,k)/u0)/u0
      Dels(4,k) = (a(3,k)*bb(4,k) - 3.0_dp*bb(3,k)/u0)*(a(1,k)*a(2,k) - 1.0_dp/u0**2) &
        & + 2.0_dp*(3.0_dp*a(1,k)*bb(2,k) - 2.0_dp*a(1,k)*bb(4,k) - 3.0_dp*bb(1,k)/u0)/u0**2
      eta(:,k) = Dels(:,k)/Del
    end do

    z1pl(:) = (eta(1,:)/2.0_dp + eta(2,:) + 5.0_dp*eta(3,:)/8.0_dp)*twopi
    z1mn(:) = (eta(1,:)/2.0_dp - eta(2,:) + 5.0_dp*eta(3,:)/8.0_dp)*twopi
    z2pl(:) = (-eta(1,:)/8.0_dp + 5.0_dp*eta(3,:)/8.0_dp + eta(4,:))*twopi
    z2mn(:) = (-eta(1,:)/8.0_dp + 5.0_dp*eta(3,:)/8.0_dp - eta(4,:))*twopi

    exptrm1(:) = exp(-min(lam1(:)*dtau(:), 35.0_dp))
    exptrm2(:) = exp(-min(lam2(:)*dtau(:), 35.0_dp))

    !! Legendre moments of the homogeneous modes
    R1(:) = -a(1,:)/lam1(:) ; R2(:) = -a(1,:)/lam2(:)
    Q1(:) = 0.5_dp*(a(1,:)*a(2,:)/(lam1(:)**2) - 1.0_dp) ; Q2(:) = 0.5_dp*(a(1,:)*a(2,:)/(lam2(:)**2) - 1.0_dp)
    S1(:) = -3.0_dp/(2.0_dp*a(4,:))*(a(1,:)*a(2,:)/lam1(:) - lam1(:))
    S2(:) = -3.0_dp/(2.0_dp*a(4,:))*(a(1,:)*a(2,:)/lam2(:) - lam2(:))

    !! Hemispheric flux (p) and second Marshak moment (q) of each mode
    p1pl(:) = (0.5_dp + R1(:) + 5.0_dp*Q1(:)/8.0_dp)*twopi
    p2pl(:) = (0.5_dp + R2(:) + 5.0_dp*Q2(:)/8.0_dp)*twopi
    q1pl(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q1(:)/8.0_dp + S1(:))*twopi
    q2pl(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q2(:)/8.0_dp + S2(:))*twopi
    p1mn(:) = (0.5_dp - R1(:) + 5.0_dp*Q1(:)/8.0_dp)*twopi
    p2mn(:) = (0.5_dp - R2(:) + 5.0_dp*Q2(:)/8.0_dp)*twopi
    q1mn(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q1(:)/8.0_dp - S1(:))*twopi
    q2mn(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q2(:)/8.0_dp - S2(:))*twopi

    f00(:) = p1mn(:)*exptrm1(:); f01(:) = p1pl(:)/exptrm1(:); f02(:) = p2mn(:)*exptrm2(:); f03(:) = p2pl(:)/exptrm2(:)
    f10(:) = q1mn(:)*exptrm1(:); f11(:) = q1pl(:)/exptrm1(:); f12(:) = q2mn(:)*exptrm2(:); f13(:) = q2pl(:)/exptrm2(:)
    f20(:) = p1pl(:)*exptrm1(:); f21(:) = p1mn(:)/exptrm1(:); f22(:) = p2pl(:)*exptrm2(:); f23(:) = p2mn(:)/exptrm2(:)
    f30(:) = q1pl(:)*exptrm1(:); f31(:) = q1mn(:)/exptrm1(:); f32(:) = q2pl(:)*exptrm2(:); f33(:) = q2mn(:)/exptrm2(:)

    expon(:) = exp(-min(tau(:)/u0, 35.0_dp))
    z1mn_up(:) = z1mn(:)*expon(2:nlev)
    z2mn_up(:) = z2mn(:)*expon(2:nlev)
    z1pl_up(:) = z1pl(:)*expon(2:nlev)
    z2pl_up(:) = z2pl(:)*expon(2:nlev)
    z1mn_down(:) = z1mn(:)*expon(1:nlay)
    z2mn_down(:) = z2mn(:)*expon(1:nlay)
    z1pl_down(:) = z1pl(:)*expon(1:nlay)
    z2pl_down(:) = z2pl(:)*expon(1:nlay)

    !! Banded matrix, entry for entry as PICASO's Mb (setM converts to LAPACK storage)
    Mb(:,:) = 0.0_dp
    B(:) = 0.0_dp

    ! Top boundary conditions
    call setM(5, 0, p1mn(1))
    call setM(5, 1, q1pl(1))
    call setM(4, 1, p1pl(1))
    call setM(4, 2, q2mn(1))
    call setM(3, 2, p2mn(1))
    call setM(3, 3, q2pl(1))
    call setM(2, 3, p2pl(1))
    call setM(6, 0, q1mn(1))

    B(1) = b_top - z1mn_down(1)
    B(2) = -b_top/4.0_dp - z2mn_down(1)

    ! Bottom boundary conditions
    call setM(5, n-2, f22(nlay))
    call setM(5, n-1, f33(nlay))
    call setM(4, n-1, f23(nlay))
    call setM(6, n-3, f21(nlay))
    call setM(6, n-2, f32(nlay))
    call setM(7, n-4, f20(nlay))
    call setM(7, n-3, f31(nlay))
    call setM(8, n-4, f30(nlay))

    B(n-1) = b_surface - z1pl_up(nlay)
    B(n) = b_surface_SH4 - z2pl_up(nlay)

    ! Fill remaining rows of matrix (layer interfaces, m is the 0-based layer above)
    do m = 0, nlay-2
      k = m + 1   ! 1-based layer above, k+1 the layer below

      call setM(5, 4*m+2, f02(k))
      call setM(5, 4*m+3, f13(k))
      call setM(5, 4*m+4, -p1pl(k+1))
      call setM(5, 4*m+5, -q1mn(k+1))

      call setM(4, 4*m+3, f03(k))
      call setM(4, 4*m+4, -q1mn(k+1))
      call setM(4, 4*m+5, -p1mn(k+1))
      call setM(4, 4*m+6, -q2pl(k+1))

      call setM(3, 4*m+4, -p1mn(k+1))
      call setM(3, 4*m+5, -q1pl(k+1))
      call setM(3, 4*m+6, -p2pl(k+1))
      call setM(3, 4*m+7, -q2mn(k+1))

      call setM(2, 4*m+5, -p1pl(k+1))
      call setM(2, 4*m+6, -q2mn(k+1))
      call setM(2, 4*m+7, -p2mn(k+1))

      call setM(1, 4*m+6, -p2mn(k+1))
      call setM(1, 4*m+7, -q2pl(k+1))

      call setM(0, 4*m+7, -p2pl(k+1))

      call setM(6, 4*m+1, f01(k))
      call setM(6, 4*m+2, f12(k))
      call setM(6, 4*m+3, f23(k))
      call setM(6, 4*m+4, -q1pl(k+1))

      call setM(7, 4*m+0, f00(k))
      call setM(7, 4*m+1, f11(k))
      call setM(7, 4*m+2, f22(k))
      call setM(7, 4*m+3, f33(k))

      call setM(8, 4*m+0, f10(k))
      call setM(8, 4*m+1, f21(k))
      call setM(8, 4*m+2, f32(k))

      call setM(9, 4*m+0, f20(k))
      call setM(9, 4*m+1, f31(k))

      call setM(10, 4*m+0, f30(k))

      B(4*m+3) = z1mn_down(k+1) - z1mn_up(k)
      B(4*m+4) = z2mn_down(k+1) - z2mn_up(k)
      B(4*m+5) = z1pl_down(k+1) - z1pl_up(k)
      B(4*m+6) = z2pl_down(k+1) - z2pl_up(k)
    end do

    call band_solve(n, Mb, B, X)

    !! Level fluxes (PICASO F X + G) - top level from layer 1, level k+1 from the bottom of layer k
    flx_down(1) = p1mn(1)*X(1) + p1pl(1)*X(2) + p2mn(1)*X(3) + p2pl(1)*X(4) + z1mn_down(1)
    flx_up(1) = p1pl(1)*X(1) + p1mn(1)*X(2) + p2pl(1)*X(3) + p2mn(1)*X(4) + z1pl_down(1)
    do k = 1, nlay
      flx_down(k+1) = f00(k)*X(4*k-3) + f01(k)*X(4*k-2) + f02(k)*X(4*k-1) + f03(k)*X(4*k) + z1mn_up(k)
      flx_up(k+1) = f20(k)*X(4*k-3) + f21(k)*X(4*k-2) + f22(k)*X(4*k-1) + f23(k)*X(4*k) + z1pl_up(k)
    end do

    !! Add the direct beam to the downward flux
    flx_down(:) = flx_down(:) + F0PI*u0*exp(-tau(:)/u0)

  contains

    !! Put PICASO's (0-based scipy) band entry (r, j) into LAPACK band storage
    subroutine setM(r, j, val)
      integer, intent(in) :: r, j
      real(dp), intent(in) :: val
      Mb(kl+1+r, j+1) = val
    end subroutine setM

  end subroutine sw_SH_four_stream

  !! Banded Gaussian elimination with partial pivoting (as LAPACK dgbsv)
  !! A(i,j) = ab(kl+ku+1+i-j,j), the top kl rows of ab are workspace
  !! (band widths kl, ku are module parameters, so the loop bounds are known at compile time)
  subroutine band_solve(n, ab, b, x)
    implicit none

    integer, intent(in) :: n
    real(dp), dimension(2*kl+ku+1, n), intent(inout) :: ab
    real(dp), dimension(n), intent(inout) :: b
    real(dp), dimension(n), intent(out) :: x

    integer :: i, j, k, p, d, jmax
    real(dp) :: amax, fac, tmp

    d = kl + ku + 1

    do k = 1, n
      jmax = min(n, k+kl+ku)

      ! Partial pivoting within the band
      p = k
      amax = abs(ab(d,k))
      do i = k+1, min(n, k+kl)
        if (abs(ab(d+i-k,k)) > amax) then
          amax = abs(ab(d+i-k,k))
          p = i
        end if
      end do
      if (p /= k) then
        do j = k, jmax
          tmp = ab(d+k-j,j)
          ab(d+k-j,j) = ab(d+p-j,j)
          ab(d+p-j,j) = tmp
        end do
        tmp = b(k) ; b(k) = b(p) ; b(p) = tmp
      end if

      ! Eliminate below the pivot
      do i = k+1, min(n, k+kl)
        fac = ab(d+i-k,k)/ab(d,k)
        if (fac == 0.0_dp) cycle
        do j = k+1, jmax
          ab(d+i-j,j) = ab(d+i-j,j) - fac*ab(d+k-j,j)
        end do
        b(i) = b(i) - fac*b(k)
      end do
    end do

    ! Back substitution
    do i = n, 1, -1
      tmp = b(i)
      do j = i+1, min(n, i+kl+ku)
        tmp = tmp - ab(d+i-j,j)*x(j)
      end do
      x(i) = tmp/ab(d,i)
    end do

  end subroutine band_solve

end module sw_SH4_mod
