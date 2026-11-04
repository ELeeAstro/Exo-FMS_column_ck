!!!
! Elspeth KH Lee - Oct 2026 : Initial version
!                - Nov 2026 : Reviewed against PICASO, provenance and differences documented;
!                             Hogan (2024) quadrature for the flux integration
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!                - Nov 2026 : Optimisations - beam integrals from the stored transmissions (one exponential
!                             per layer and angle), compile-time band widths in the solver,
!                             reciprocals, contiguous per g-ordinance inputs
!
! lw: Four-stream spherical harmonic method (SH4) with multiple scattering
!     (Rooney et al. 2024), with the level up and down fluxes found using the source
!     function technique (Toon et al. 1989)
!     Taken from PICASO (Batalha et al. 2019; Rooney et al. 2024) and converted to Fortran:
!     fluxes.get_thermal_SH with stream = 4 and fluxes.setup_4_stream_fluxes (four-stream system),
!     checked line by line against PICASO commit 4524b70 (Oct 2026). The banded matrix, the
!     linear-in-tau thermal particular solution, the delta-M scaling (f = g**4, as PICASO's
!     optics) and the source-function intensity integration follow PICASO term by term.
!!!

module lw_SH4_mod
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
  integer, parameter :: nstr = 4

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

  !! Band widths of the four-stream matrix
  integer, parameter :: kl = 5, ku = 5

  private :: lw_SH_four_stream, band_solve, e_alp, e_bet
  public :: lw_SH4

contains

  subroutine lw_SH4(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
        call lw_SH_four_stream(nlay, nlev, be_b, be_int(b), tau_g, ssa_g, gg_g, lw_up_g(:), lw_down_g(:))
        lw_up(:) = lw_up(:) + lw_up_g(:) * gw(g)
        lw_down(:) = lw_down(:) + lw_down_g(:) * gw(g)
      end do
    end do

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_SH4

  subroutine lw_SH_four_stream(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables
    integer :: k, m, n, r, c, row, col0, i0
    real(dp) :: tautop, Btop, u, ru, P2u, P3u
    real(dp), dimension(nlay) :: dtau, w0, ff, B0, B1
    real(dp), dimension(nlay) :: w_multi0, w_multi1, w_multi2, w_multi3, a0, a1, a2, a3
    real(dp), dimension(nlay) :: beta, gama, disc, lam1, lam2, e1, e2
    real(dp), dimension(nlay) :: R1, R2, Q1, Q2, S1, S2
    real(dp), dimension(nlay) :: p1pl, p2pl, q1pl, q2pl, p1mn, p2mn, q1mn, q2mn
    real(dp), dimension(4,4,nlay) :: Ftop, Fbot
    real(dp), dimension(4,nlay) :: zdown, zup
    real(dp), dimension(2) :: bsurf
    real(dp), dimension(nlay) :: X1, X2, X3, X4
    real(dp), dimension(nlay) :: ex, ev1, od1, ev2, od2, T_up, T_dn, Ip1, Nup, Ndn, Aup, Adn
    real(dp), dimension(nlay) :: ea1, eb1, ea2, eb2
    real(dp), dimension(nlay) :: src_up, src_dn, el1, el2, wa0
    real(dp), dimension(4*nlay) :: rhs, X
    real(dp), dimension(2*kl+ku+1, 4*nlay) :: Am
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    n = 4*nlay

    !! Delta-M scaling (PICASO convention, f = g**nstr)
    ff(:) = g_in(:)**nstr
    w0(:) = w_in(:)*(1.0_dp - ff(:))/(1.0_dp - w_in(:)*ff(:))
    w0(:) = min(w0(:), 1.0_dp - 1.0e-6_dp)
    dtau(:) = (tau_in(2:nlev) - tau_in(1:nlay))*(1.0_dp - w_in(:)*ff(:))

    !! Phase function Legendre coefficients and the PICASO 'a' coefficients
    w_multi0(:) = 1.0_dp
    w_multi1(:) = 3.0_dp*(g_in(:) - ff(:))/(1.0_dp - ff(:))
    w_multi2(:) = 5.0_dp*(g_in(:)**2 - ff(:))/(1.0_dp - ff(:))
    w_multi3(:) = 7.0_dp*(g_in(:)**3 - ff(:))/(1.0_dp - ff(:))
    a0(:) = 1.0_dp - w0(:)*w_multi0(:)
    a1(:) = 3.0_dp - w0(:)*w_multi1(:)
    a2(:) = 5.0_dp - w0(:)*w_multi2(:)
    a3(:) = 7.0_dp - w0(:)*w_multi3(:)

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

    !! Eigenvalues; lam2 in the cancellation free form lam1**2 * lam2**2 = gama
    beta(:) = a0(:)*a1(:) + 4.0_dp*a0(:)*a3(:)/9.0_dp + a2(:)*a3(:)/9.0_dp
    gama(:) = a0(:)*a1(:)*a2(:)*a3(:)/9.0_dp
    disc(:) = sqrt(max(beta(:)**2 - 4.0_dp*gama(:), 0.0_dp))
    lam1(:) = sqrt((beta(:) + disc(:))/2.0_dp)
    lam2(:) = sqrt(2.0_dp*gama(:)/(beta(:) + disc(:)))

    e1(:) = exp(-min(lam1(:)*dtau(:), 35.0_dp))
    e2(:) = exp(-min(lam2(:)*dtau(:), 35.0_dp))
    ! uncapped exp(-lam dtau) for the beam integrals (= e1, e2 unless lam dtau > 35)
    el1(:) = e1(:)
    where (lam1(:)*dtau(:) > 35.0_dp) el1(:) = exp(-lam1(:)*dtau(:))
    el2(:) = e2(:)
    where (lam2(:)*dtau(:) > 35.0_dp) el2(:) = exp(-lam2(:)*dtau(:))

    !! Legendre moments of the homogeneous modes
    R1(:) = -a0(:)/lam1(:) ; R2(:) = -a0(:)/lam2(:)
    Q1(:) = 0.5_dp*(a0(:)*a1(:)/lam1(:)**2 - 1.0_dp) ; Q2(:) = 0.5_dp*(a0(:)*a1(:)/lam2(:)**2 - 1.0_dp)
    S1(:) = -3.0_dp/(2.0_dp*a3(:))*(a0(:)*a1(:)/lam1(:) - lam1(:))
    S2(:) = -3.0_dp/(2.0_dp*a3(:))*(a0(:)*a1(:)/lam2(:) - lam2(:))

    !! Hemispheric flux (p) and second Marshak moment (q) of each mode
    p1pl(:) = twopi*(0.5_dp + R1(:) + 5.0_dp*Q1(:)/8.0_dp)
    p2pl(:) = twopi*(0.5_dp + R2(:) + 5.0_dp*Q2(:)/8.0_dp)
    q1pl(:) = twopi*(-1.0_dp/8.0_dp + 5.0_dp*Q1(:)/8.0_dp + S1(:))
    q2pl(:) = twopi*(-1.0_dp/8.0_dp + 5.0_dp*Q2(:)/8.0_dp + S2(:))
    p1mn(:) = twopi*(0.5_dp - R1(:) + 5.0_dp*Q1(:)/8.0_dp)
    p2mn(:) = twopi*(0.5_dp - R2(:) + 5.0_dp*Q2(:)/8.0_dp)
    q1mn(:) = twopi*(-1.0_dp/8.0_dp + 5.0_dp*Q1(:)/8.0_dp - S1(:))
    q2mn(:) = twopi*(-1.0_dp/8.0_dp + 5.0_dp*Q2(:)/8.0_dp - S2(:))

    !! Moments at the layer top: rows (down 1st, down 2nd, up 1st, up 2nd),
    !! columns modes (e^-lam1 t, e^+lam1 t, e^-lam2 t, e^+lam2 t)
    do k = 1, nlay
      Ftop(1,:,k) = (/p1mn(k), p1pl(k), p2mn(k), p2pl(k)/)
      Ftop(2,:,k) = (/q1mn(k), q1pl(k), q2mn(k), q2pl(k)/)
      Ftop(3,:,k) = (/p1pl(k), p1mn(k), p2pl(k), p2mn(k)/)
      Ftop(4,:,k) = (/q1pl(k), q1mn(k), q2pl(k), q2mn(k)/)
      Fbot(:,1,k) = Ftop(:,1,k)*e1(k)
      Fbot(:,2,k) = Ftop(:,2,k)/e1(k)
      Fbot(:,3,k) = Ftop(:,3,k)*e2(k)
      Fbot(:,4,k) = Ftop(:,4,k)/e2(k)
    end do

    !! Particular solution moments at layer top (down) and bottom (up), linear in tau Planck
    zdown(1,:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp - B1(:)/a1(:))
    zdown(2,:) = -twopi*0.5_dp*(1.0_dp - w0(:))/(4.0_dp*a0(:)) * B0(:)
    zdown(3,:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp + B1(:)/a1(:))
    zdown(4,:) = zdown(2,:)
    zup(1,:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp - B1(:)/a1(:) + B1(:)*dtau(:)/2.0_dp)
    zup(2,:) = -twopi*0.5_dp*(1.0_dp - w0(:))/(4.0_dp*a0(:)) * (B0(:) + B1(:)*dtau(:))
    zup(3,:) = twopi*(1.0_dp - w0(:))/a0(:) * (B0(:)/2.0_dp + B1(:)/a1(:) + B1(:)*dtau(:)/2.0_dp)
    zup(4,:) = zup(2,:)

    !! Top boundary, isothermal layer above the model top (as lw_Toon89)
    tautop = dtau(1)*exp(-1.0_dp)
    Btop = (1.0_dp - exp(-tautop/ubari))*be(1)

    !! Build the banded matrix, LAPACK storage A(i,j) = Am(kl+ku+1+i-j,j)
    !! Unknowns per layer k: X(4k-3:4k)
    Am(:,:) = 0.0_dp
    i0 = kl + ku + 1

    ! Top boundary: downward moments = those of the incident isotropic intensity
    do r = 1, 2
      do c = 1, 4
        Am(i0+r-c,c) = Ftop(r,c,1)
      end do
    end do
    rhs(1) = pi*Btop - zdown(1,1)
    rhs(2) = -pi*Btop/4.0_dp - zdown(2,1)

    ! Layer interfaces: continuity of all four moments
    do k = 1, nlay-1
      col0 = 4*(k-1)
      do r = 1, 4
        row = 2 + col0 + r
        do c = 1, 4
          Am(i0+row-(col0+c),col0+c) = Fbot(r,c,k)
          Am(i0+row-(col0+4+c),col0+4+c) = -Ftop(r,c,k+1)
        end do
        rhs(row) = zdown(r,k+1) - zup(r,k)
      end do
    end do

    ! Bottom boundary: gas giant internal flux, I_up = I_down + B_int (unit return)
    bsurf(1) = pi*be_int
    bsurf(2) = -pi*be_int/4.0_dp
    col0 = n - 4
    do r = 1, 2
      row = n - 2 + r
      do c = 1, 4
        Am(i0+row-(col0+c),col0+c) = Fbot(2+r,c,nlay) - Fbot(r,c,nlay)
      end do
      rhs(row) = bsurf(r) - zup(2+r,nlay) + zup(r,nlay)
    end do

    call band_solve(n, Am, rhs, X)

    !! Growing mode coefficients rescaled to their value at the layer bottom (avoids overflow)
    X1(:) = X(1:n:4)
    X2(:) = X(2:n:4)/e1(:)
    X3(:) = X(3:n:4)
    X4(:) = X(4:n:4)/e2(:)

    !! Particular solution first moment (constant in tau)
    Ip1(:) = (1.0_dp - w0(:))/(a0(:)*a1(:))*B1(:)
    wa0(:) = w_multi0(:)*(1.0_dp - w0(:))/a0(:)

    ! Zero the total flux arrays
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp

    do m = 1, nmu

      u = uarr(m)
      ru = ruarr(m)
      P2u = (3.0_dp*u**2 - 1.0_dp)/2.0_dp
      P3u = (5.0_dp*u**3 - 3.0_dp*u)/2.0_dp
      ex(:) = exp(-dtau(:)*ru)

      !! Scattering phase coefficients split into even and odd Legendre parts,
      !! at -u (and for the growing mode of each pair) the odd part flips sign
      ev1(:) = w_multi0(:) + w_multi2(:)*P2u*Q1(:)
      od1(:) = w_multi1(:)*u*R1(:) + w_multi3(:)*P3u*S1(:)
      ev2(:) = w_multi0(:) + w_multi2(:)*P2u*Q2(:)
      od2(:) = w_multi1(:)*u*R2(:) + w_multi3(:)*P3u*S2(:)

      ea1(:) = e_alp(nlay, ru + lam1(:), ex(:), el1(:))
      eb1(:) = e_bet(nlay, lam1(:), ru, dtau(:), ex(:), el1(:))
      ea2(:) = e_alp(nlay, ru + lam2(:), ex(:), el2(:))
      eb2(:) = e_bet(nlay, lam2(:), ru, dtau(:), ex(:), el2(:))

      !! Thermal particular solution integrated along the beam
      T_up(:) = u*(B0(:)*(1.0_dp - ex(:)) + B1(:)*(u - (dtau(:) + u)*ex(:)))
      T_dn(:) = u*(B0(:)*(1.0_dp - ex(:)) + B1(:)*(dtau(:) - u + u*ex(:)))
      Nup(:) = wa0(:)*T_up(:) + w_multi1(:)*u*Ip1(:)*u*(1.0_dp - ex(:))
      Ndn(:) = wa0(:)*T_dn(:) - w_multi1(:)*u*Ip1(:)*u*(1.0_dp - ex(:))

      !! Homogeneous solution integrated along the beam
      Aup(:) = X1(:)*(ev1(:) + od1(:))*ea1(:) + X2(:)*(ev1(:) - od1(:))*eb1(:) &
        & + X3(:)*(ev2(:) + od2(:))*ea2(:) + X4(:)*(ev2(:) - od2(:))*eb2(:)
      Adn(:) = X1(:)*(ev1(:) - od1(:))*eb1(:) + X2(:)*(ev1(:) + od1(:))*ea1(:) &
        & + X3(:)*(ev2(:) - od2(:))*eb2(:) + X4(:)*(ev2(:) + od2(:))*ea2(:)

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

  end subroutine lw_SH_four_stream

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

end module lw_SH4_mod
