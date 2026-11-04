!!!
! Elspeth KH Lee - Oct 2026 : Initial version
!                - Nov 2026 : Reviewed against PICASO, provenance and differences documented;
!                             Hogan (2024) quadrature for the flux integration
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!                - Nov 2026 : Optimisations - beam integrals from the stored transmissions (one exponential
!                             per layer and angle), compile-time band widths in the solver,
!                             reciprocals, contiguous per g-ordinance inputs
!
! lw: Two-stream spherical harmonic method (SH2) with multiple scattering
!     (Rooney et al. 2024), with the level up and down fluxes found using the source
!     function technique (Toon et al. 1989)
!     Taken from PICASO (Batalha et al. 2019; Rooney et al. 2024) and converted to Fortran:
!     fluxes.get_thermal_SH with stream = 2 and fluxes.setup_2_stream_fluxes (two-stream, tridiagonal-pair system),
!     checked line by line against PICASO commit 4524b70 (Oct 2026). The banded matrix, the
!     linear-in-tau thermal particular solution, the delta-M scaling (f = g**2, as PICASO's
!     optics) and the source-function intensity integration follow PICASO term by term.
!!!

module lw_SH2_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4
  use BB_integrate_mod, only : BB_integrate
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Number of streams and the Legendre order of the delta-M scaling
  integer, parameter :: nstr = 2

  real(dp), parameter :: ubari = 0.5_dp

  !! Gauss quadrature variables, cosine angle values (uarr) and weights (wuarr)
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these

  !! Gauss–Jacobi-5 quadrature for 4 nodes - 8 stream (Hogan 2024)
  !! Hogan's weights w are for the measure 2 mu dmu (sum 1, F = pi sum w I); here the 2 pi
  !! is carried in the intensities (F = 2 pi sum wuarr I), so wuarr = w/2 (sum 0.5)
  integer, parameter :: nmu = 4
  real(dp), dimension(nmu), parameter :: uarr = (/0.0454586727_dp, 0.2322334416_dp, 0.5740198775_dp, 0.9030775973_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.0092068785_dp, 0.1285704278_dp, 0.4323381850_dp, 0.4298845087_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = 0.5_dp * w
  real(dp), dimension(nmu), parameter :: ruarr = 1.0_dp/uarr

  !! Band widths of the two-stream matrix
  integer, parameter :: kl = 2, ku = 2

  private :: lw_SH_two_stream, band_solve, e_alp, e_bet
  public :: lw_SH2

contains

  subroutine lw_SH2(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
    real(dp), dimension(nlev) :: be_b, tau_g
    real(dp), dimension(nlay) :: ssa_g, gg_g

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
      be_b(:) = be(b,:)
      do g = 1, ng
        ! Contiguous copies (avoids array temporaries for the strided slices)
        tau_g(:) = tau_e(g,b,:)
        ssa_g(:) = ssa(g,b,:)
        gg_g(:) = gg(g,b,:)
        call lw_SH_two_stream(nlay, nlev, be_b, be_int(b), tau_g, ssa_g, gg_g, lw_up_g(:), lw_down_g(:))
        lw_up(:) = lw_up(:) + lw_up_g(:) * gw(g)
        lw_down(:) = lw_down(:) + lw_down_g(:) * gw(g)
      end do
    end do

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_SH2

  subroutine lw_SH_two_stream(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables
    integer :: k, m, n, i0
    real(dp) :: tautop, Btop, u, ru
    real(dp), dimension(nlay) :: dtau, w0, ff, B0, B1
    real(dp), dimension(nlay) :: w_multi0, w_multi1, a0, a1
    real(dp), dimension(nlay) :: lam, q, Q1, Q2, exptrm
    real(dp), dimension(nlay) :: zmn_down, zmn_up, zpl_down, zpl_up
    real(dp), dimension(nlay) :: X1, X2
    real(dp), dimension(nlay) :: ex, c_ev, c_od, T_up, T_dn, Ip1, Nup, Ndn, Aup, Adn, ea, eb
    real(dp), dimension(nlay) :: src_up, src_dn, el, wa0
    real(dp), dimension(2*nlay) :: rhs, X
    real(dp), dimension(2*kl+ku+1, 2*nlay) :: Am
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    n = 2*nlay

    !! Delta-M scaling (PICASO convention, f = g**nstr)
    ff(:) = g_in(:)**nstr
    w0(:) = w_in(:)*(1.0_dp - ff(:))/(1.0_dp - w_in(:)*ff(:))
    w0(:) = min(w0(:), 1.0_dp - 1.0e-6_dp)
    dtau(:) = (tau_in(2:nlev) - tau_in(1:nlay))*(1.0_dp - w_in(:)*ff(:))

    !! Phase function Legendre coefficients and the PICASO 'a' coefficients
    w_multi0(:) = 1.0_dp
    w_multi1(:) = 3.0_dp*(g_in(:) - ff(:))/(1.0_dp - ff(:))
    a0(:) = 1.0_dp - w0(:)*w_multi0(:)
    a1(:) = 3.0_dp - w0(:)*w_multi1(:)

    do k = 1, nlay
      if (dtau(k) <= 1.0e-6_dp) then
        ! For low optical depths use the isothermal approimation
        B1(k) = 0.0_dp
        B0(k) = 0.5_dp*(be(k+1) + be(k))
      else
        B1(k) = (be(k+1) - be(k))/dtau(k) ! Linear in tau term
        B0(k) = be(k)
      endif
    end do

    !! Homogeneous solution: fluxes of mode 1 (decays with depth) and mode 2 (grows)
    lam(:) = sqrt(a0(:)*a1(:))
    q(:) = lam(:)/a1(:)
    Q1(:) = twopi*(0.5_dp + q(:))
    Q2(:) = twopi*(0.5_dp - q(:))
    exptrm(:) = exp(-min(lam(:)*dtau(:), 35.0_dp))
    ! uncapped exp(-lam dtau) for the beam integrals (= exptrm unless lam dtau > 35)
    el(:) = exptrm(:)
    where (lam(:)*dtau(:) > 35.0_dp) el(:) = exp(-lam(:)*dtau(:))

    !! Particular solution fluxes at layer top (down) and bottom (up), linear in tau Planck
    zmn_down(:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp - B1(:)/a1(:))
    zmn_up(:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp - B1(:)/a1(:) + B1(:)*dtau(:)/2.0_dp)
    zpl_down(:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp + B1(:)/a1(:))
    zpl_up(:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp + B1(:)/a1(:) + B1(:)*dtau(:)/2.0_dp)

    !! Top boundary, isothermal layer above the model top (as lw_Toon89)
    tautop = dtau(1)*exp(-1.0_dp)
    Btop = (1.0_dp - exp(-tautop/ubari))*be(1)

    !! Build the banded matrix, LAPACK storage A(i,j) = Am(kl+ku+1+i-j,j)
    !! Unknowns per layer k: X(2k-1) (mode 1), X(2k) (mode 2)
    Am(:,:) = 0.0_dp
    i0 = kl + ku + 1

    ! Top boundary: downward flux = incident isotropic flux
    Am(i0+1-1,1) = Q1(1) ; Am(i0+1-2,2) = Q2(1)
    rhs(1) = pi*Btop - zmn_down(1)

    ! Layer interfaces: continuity of the downward and upward fluxes
    do k = 1, nlay-1
      m = 2*k
      ! downward flux
      Am(i0+m-(m-1),m-1) = Q1(k)*exptrm(k)
      Am(i0+m-m,m) = Q2(k)/exptrm(k)
      Am(i0+m-(m+1),m+1) = -Q1(k+1)
      Am(i0+m-(m+2),m+2) = -Q2(k+1)
      rhs(m) = zmn_down(k+1) - zmn_up(k)
      ! upward flux
      Am(i0+(m+1)-(m-1),m-1) = Q2(k)*exptrm(k)
      Am(i0+(m+1)-m,m) = Q1(k)/exptrm(k)
      Am(i0+(m+1)-(m+1),m+1) = -Q2(k+1)
      Am(i0+(m+1)-(m+2),m+2) = -Q1(k+1)
      rhs(m+1) = zpl_down(k+1) - zpl_up(k)
    end do

    ! Bottom boundary: gas giant internal flux, I_up = I_down + B_int (unit return)
    Am(i0+n-(n-1),n-1) = Q2(nlay)*exptrm(nlay) - Q1(nlay)*exptrm(nlay)
    Am(i0+n-n,n) = Q1(nlay)/exptrm(nlay) - Q2(nlay)/exptrm(nlay)
    rhs(n) = pi*be_int - zpl_up(nlay) + zmn_up(nlay)

    call band_solve(n, Am, rhs, X)

    !! Mode 2 coefficient rescaled to its value at the layer bottom (avoids overflow)
    X1(:) = X(1:n:2)
    X2(:) = X(2:n:2)/exptrm(:)

    !! Particular solution first moment (constant in tau)
    Ip1(:) = (1.0_dp - w0(:))/(a0(:)*a1(:))*B1(:)
    wa0(:) = w_multi0(:)*(1.0_dp - w0(:))/a0(:)

    ! Zero the total flux arrays
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp

    do m = 1, nmu

      u = uarr(m)
      ru = ruarr(m)
      ex(:) = exp(-dtau(:)*ru)

      !! Scattering phase coefficients for each mode at +u; at -u the odd terms flip
      !! mode 1 Legendre moments (1, -q), mode 2 (1, q)
      c_ev(:) = w_multi0(:)
      c_od(:) = -w_multi1(:)*u*q(:)

      ea(:) = e_alp(nlay, ru + lam(:), ex(:), el(:))
      eb(:) = e_bet(nlay, lam(:), ru, dtau(:), ex(:), el(:))

      !! Thermal particular solution integrated along the beam
      T_up(:) = u*(B0(:)*(1.0_dp - ex(:)) + B1(:)*(u - (dtau(:) + u)*ex(:)))
      T_dn(:) = u*(B0(:)*(1.0_dp - ex(:)) + B1(:)*(dtau(:) - u + u*ex(:)))
      Nup(:) = wa0(:)*T_up(:) + w_multi1(:)*u*Ip1(:)*u*(1.0_dp - ex(:))
      Ndn(:) = wa0(:)*T_dn(:) - w_multi1(:)*u*Ip1(:)*u*(1.0_dp - ex(:))

      !! Homogeneous solution integrated along the beam
      Aup(:) = X1(:)*(c_ev(:) + c_od(:))*ea(:) + X2(:)*(c_ev(:) - c_od(:))*eb(:)
      Adn(:) = X1(:)*(c_ev(:) - c_od(:))*eb(:) + X2(:)*(c_ev(:) + c_od(:))*ea(:)

      src_up(:) = twopi*ru*(w0(:)*(Aup(:) + Nup(:)) + (1.0_dp - w0(:))*T_up(:))
      src_dn(:) = twopi*ru*(w0(:)*(Adn(:) + Ndn(:)) + (1.0_dp - w0(:))*T_dn(:))

      !! Perform downward loop first
      ! Top boundary condition - intensity downward from top boundary (tautop, assumed isothermal)
      lw_down_g(1) = twopi*(1.0_dp - exp(-tautop*ru))*be(1)
      do k = 1, nlay
        lw_down_g(k+1) = lw_down_g(k)*ex(k) + src_dn(k)
      end do

      lw_up_g(nlev) = lw_down_g(nlev) + twopi*be_int
      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*ex(k) + src_up(k)
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      flx_down(:) = flx_down(:) + lw_down_g(:) * wuarr(m)
      flx_up(:) = flx_up(:) + lw_up_g(:) * wuarr(m)

    end do

  end subroutine lw_SH_two_stream

  !! int_0^dtau exp(-alp t) dt, alp = 1/u + lam > 0, from ex = exp(-dtau/u) and el = exp(-lam dtau)
  function e_alp(nlay, alp, ex, el) result(res)
    implicit none
    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: alp, ex, el
    real(dp), dimension(nlay) :: res

    res(:) = (1.0_dp - ex(:)*el(:))/alp(:)

  end function e_alp

  !! (exp(-lam dtau) - exp(-dtau/u))/(1/u - lam), with the series form as 1/u - lam -> 0,
  !! from ex = exp(-dtau/u) and el = exp(-lam dtau)
  function e_bet(nlay, lam, ru, dtau, ex, el) result(res)
    implicit none
    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: lam, dtau, ex, el
    real(dp), intent(in) :: ru
    real(dp), dimension(nlay) :: res

    integer :: k
    real(dp) :: bet, x

    do k = 1, nlay
      bet = ru - lam(k)
      x = bet*dtau(k)
      if (abs(x) < 1.0e-6_dp) then
        res(k) = dtau(k)*ex(k)*(1.0_dp + x/2.0_dp + x*x/6.0_dp)
      else
        res(k) = (el(k) - ex(k))/bet
      end if
    end do

  end function e_bet

  !! Banded Gaussian elimination with partial pivoting (as LAPACK dgbsv)
  !! A(i,j) = ab(kl+ku+1+i-j,j), the top kl rows of ab are fill-in workspace
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

end module lw_SH2_mod
