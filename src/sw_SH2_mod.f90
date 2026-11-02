!!!
! Elspeth KH Lee - Jun 2024 : Initial version
!                - Nov 2026 : Rewritten following PICASO exactly, enabled
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!
! sw: Two-stream spherical harmonic method (SH2) with multiple scattering (Rooney et al. 2024)
!     Taken from PICASO (Batalha et al. 2019; Rooney et al. 2024) and converted to Fortran:
!     fluxes.get_reflected_SH and fluxes.setup_2_stream_fluxes, checked against PICASO commit 4524b70
!     (Oct 2026) to 1e-13. Delta-M scaling (f = g**2) and one-term HG phase function as PICASO.
!     Differences to PICASO: no Rayleigh, gas giant lower boundary, own banded solver,
!     fluxes found row by row, single zenith angle (no zcorr). None change the solution.
!     TTHG option: as PICASO, but delta-M uses the TTHG chi_nstr for w0 and dtau too (consistent).
!     Pros: Fast multiple scattering method (gives the same answer as sw_adding)
!     Cons: Two-stream only, slower than sw_adding
!!!

module sw_SH2_mod
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
  integer, parameter :: nstr = 2

  !! Use Two-Term HG function (TTHG), parameters are the PICASO defaults
  !! forward fraction = frac_a + frac_b*g_back**frac_c, g_forward = c_forward*g, g_back = c_back*g
  logical, parameter :: TTHG = .False.
  real(dp), parameter :: tthg_forward = 1.0_dp, tthg_back = -0.5_dp
  real(dp), parameter :: tthg_frac_a = 1.0_dp, tthg_frac_b = -1.0_dp
  integer, parameter :: tthg_frac_c = 2

  private :: sw_SH_two_stream, band_solve
  public :: sw_SH2

contains

  subroutine sw_SH2(nlay, nlev, nb, ng, gw, tau_e, mu_z, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
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
          call sw_SH_two_stream(nlay, nlev, Finc(b), mu_z(nlev), tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), &
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

  end subroutine sw_SH2

  subroutine sw_SH_two_stream(nlay, nlev, F0PI, u0, tau_in, w_in, g_in, flx_down, flx_up)
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
    integer, parameter :: kl = 2, ku = 2
    real(dp) :: b_top, b_surface
    real(dp), dimension(nstr) :: Pu0
    real(dp), dimension(nlay) :: dtau, w0, f_deltaM, g_f, g_b, f_f
    real(dp), dimension(nstr,nlay) :: chi
    real(dp), dimension(nstr,nlay) :: w_multi, a, bb, eta
    real(dp), dimension(nlay) :: Del, lam, exptrm, q, Q1, Q2, Q1mn, Q2mn, Q1pl, Q2pl
    real(dp), dimension(nlay) :: zmn, zpl, zmn_up, zpl_up, zmn_down, zpl_down
    real(dp), dimension(nlev) :: tau, expon
    real(dp), dimension(2*kl+ku+1, 2*nlay) :: Mb
    real(dp), dimension(2*nlay) :: B, X

    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w_in(:) <= 1.0e-12_dp)) then
      flx_down(:) = F0PI * u0 * exp(-tau_in(:)/u0)
      flx_up(:) = 0.0_dp
      return
    end if

    n = 2*nlay

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

    !! Direct beam particular solution
    Del(:) = (1.0_dp/u0)**2 - a(1,:)*a(2,:)
    eta(1,:) = (bb(2,:)/u0 - a(2,:)*bb(1,:))/Del(:)
    eta(2,:) = (bb(1,:)/u0 - a(1,:)*bb(2,:))/Del(:)

    lam(:) = sqrt(a(1,:)*a(2,:))
    exptrm(:) = exp(-min(lam(:)*dtau(:), 35.0_dp))

    !! Parameters in the matrices
    q(:) = lam(:)/a(2,:)
    Q1(:) = (0.5_dp + q(:))*twopi
    Q2(:) = (0.5_dp - q(:))*twopi

    Q1mn(:) = Q1(:)*exptrm(:);  Q2mn(:) = Q2(:)*exptrm(:)
    Q1pl(:) = Q1(:)/exptrm(:);  Q2pl(:) = Q2(:)/exptrm(:)

    zmn(:) = (0.5_dp*eta(1,:) - eta(2,:))*twopi
    zpl(:) = (0.5_dp*eta(1,:) + eta(2,:))*twopi
    expon(:) = exp(-tau(:)/u0)
    zmn_up(:) = zmn(:)*expon(2:nlev)
    zpl_up(:) = zpl(:)*expon(2:nlev)
    zmn_down(:) = zmn(:)*expon(1:nlay)
    zpl_down(:) = zpl(:)*expon(1:nlay)

    !! Banded matrix, entry for entry as PICASO's Mb (setM converts to LAPACK storage)
    Mb(:,:) = 0.0_dp
    B(:) = 0.0_dp

    ! First row: BC 1
    call setM(2, 0, Q1(1))
    call setM(1, 1, Q2(1))
    B(1) = b_top - zmn_down(1)

    ! Last row: BC 4
    call setM(3, n-2, Q2mn(nlay))
    call setM(2, n-1, Q1pl(nlay))
    B(n) = b_surface - zpl_up(nlay)

    ! Remaining rows
    do m = 1, nlay-1
      ! Layer below the interface
      call setM(0, 2*m+1, -Q2(m+1))
      call setM(1, 2*m, -Q1(m+1))
      call setM(1, 2*m+1, -Q1(m+1))
      call setM(2, 2*m, -Q2(m+1))
      ! Layer above the interface
      call setM(2, 2*m-1, Q2pl(m))
      call setM(3, 2*m-2, Q1mn(m))
      call setM(3, 2*m-1, Q1pl(m))
      call setM(4, 2*m-2, Q2mn(m))
      B(2*m) = zmn_down(m+1) - zmn_up(m)
      B(2*m+1) = zpl_down(m+1) - zpl_up(m)
    end do

    call band_solve(n, kl, ku, Mb, B, X)

    !! Level fluxes (PICASO F X + G) - top level from layer 1, level k+1 from the bottom of layer k
    flx_down(1) = Q1(1)*X(1) + Q2(1)*X(2) + zmn_down(1)
    flx_up(1) = Q2(1)*X(1) + Q1(1)*X(2) + zpl_down(1)
    do k = 1, nlay
      flx_down(k+1) = Q1mn(k)*X(2*k-1) + Q2pl(k)*X(2*k) + zmn_up(k)
      flx_up(k+1) = Q2mn(k)*X(2*k-1) + Q1pl(k)*X(2*k) + zpl_up(k)
    end do

    !! Add the direct beam to the downward flux
    flx_down(:) = flx_down(:) + F0PI*u0*expon(:)

  contains

    !! Put PICASO's (0-based scipy) band entry (r, j) into LAPACK band storage
    subroutine setM(r, j, val)
      integer, intent(in) :: r, j
      real(dp), intent(in) :: val
      Mb(kl+1+r, j+1) = val
    end subroutine setM

  end subroutine sw_SH_two_stream

  !! Banded Gaussian elimination with partial pivoting (as LAPACK dgbsv)
  !! A(i,j) = ab(kl+ku+1+i-j,j), the top kl rows of ab are workspace
  subroutine band_solve(n, kl, ku, ab, b, x)
    implicit none

    integer, intent(in) :: n, kl, ku
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

end module sw_SH2_mod
