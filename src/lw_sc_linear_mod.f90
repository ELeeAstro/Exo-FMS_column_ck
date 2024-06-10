!!!
! Elspeth KH Lee - May 2021 : Initial version
!                - Oct 2021 : adding method & Bezier interpolation
!                - Aug 2023 : Change quadrature following Li (2000)
!                - Jun 2024 : Change quadrature following Hogan (2004)
! lw: Two-stream method following the short characteristics method (e.g. Helios-r2: Kitzmann et al. 2018)
!     Uses the method of short characteristics (Olson & Kunasz 1987) with linear interpolants.
!     Pros: Very fast, accurate at high optical depths, very stable
!     Cons: No lw scattering (but could combine with AA if needed)
!!!

module lw_sc_linear_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4  
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp
  real(dp), parameter :: hp = 6.62607015e-34_dp
  real(dp), parameter :: kb = 1.380649e-23_dp
  real(dp), parameter :: c_s = 2.99792458e8_dp
  real(dp), parameter :: c1 = (hp * c_s) / kb
  real(dp), parameter :: c2 = c_s**2
  real(dp), parameter :: n2 = 2.0_dp * hp * c2

  !! Gauss quadrature variables, cosine angle values (uarr) and weights (w)
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these
  
  !! Optimised quadrature for 1 node (Hogan 2024)
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.6096748751_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)

  !! Gauss–Jacobi-5 quadrature for 2 nodes (Hogan 2024)
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.2509907356_dp, 0.7908473988_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.2300253764_dp, 0.7699746236_dp/)

  private :: lw_shortchar_linear, BB_integrate
  public :: lw_sc_linear

contains


  subroutine lw_sc_linear(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, Tint, lw_up, lw_down, lw_net, olr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev, nb, ng
    real(dp), dimension(ng), intent(in) :: gw
    real(dp), dimension(nb+1), intent(in) :: wn_e
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(ng,nb,nlev), intent(in) :: tau_e
    real(dp), intent(in) :: Tint

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down, lw_net
    real(dp), intent(out) :: olr

    !! Work variables
    integer :: k, b, g
    real(dp), dimension(nlev) :: Te
    real(dp), dimension(nb,nlev) :: be
    real(dp), dimension(nb) :: be_int
    real(dp), dimension(nb,nlev) :: lw_down_b, lw_up_b
    real(dp), dimension(ng,nlev) :: lw_down_g, lw_up_g

    !! Use WENO4 method to (smoothly) interpolate layers to levels
    Te(:) = interpolate_weno4(pe, pl, Tl, .False.)

    !! Edges are linearly interpolated to avoid overshoot
    Te(1) = 10.0_dp**(log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Tl(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    !! Find integrated planck function for each band for each level
    do k = 1, nlev
      call BB_integrate(nb, Te(k), wn_e, be(:,k))   ! Integrated planck function intensity at levels for each band
    end do
    call BB_integrate(nb, Tint, wn_e, be_int) ! Integrated planck function intensity for internal temperature

    !! Longwave flux calculation
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp
    do b = 1, nb
      lw_up_b(b,:) = 0.0_dp
      lw_down_b(b,:) = 0.0_dp
      do g = 1, ng
        call lw_shortchar_linear(nlay, nlev, be(b,:), be_int(b), tau_e(g,b,:), lw_up_g(g,:), lw_down_g(g,:))
        lw_up_b(b,:) = lw_up_b(b,:) + lw_up_g(g,:) * gw(g)
        lw_down_b(b,:) = lw_down_b(b,:) + lw_down_g(g,:) * gw(g)
      end do
      lw_up(:) = lw_up(:) + lw_up_b(b,:)
      lw_down(:) = lw_down(:) + lw_down_b(b,:)
    end do

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_sc_linear

  subroutine lw_shortchar_linear(nlay, nlev, be, be_int, tau_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables and arrays
    integer :: k, m
    real(dp), dimension(nlay) :: dtau, edel
    real(dp), dimension(nlay) :: del, e0i, e1i, e1i_del
    real(dp), dimension(nlay) :: Am, Bm, Gp, Bp
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    !! Calculate dtau in each layer
    dtau(:) = tau_in(2:) - tau_in(1:nlay)

    ! Zero the total flux arrays
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nmu

      del(:) = dtau(:)/uarr(m)
      edel(:) = exp(-del(:))
      e0i(:) = 1.0_dp - edel(:)

      !! Prepare loop
      ! Olson & Kunasz (1987) linear interpolant parameters
      where (edel(:) > 0.999_dp)
        ! If we are in very low optical depth regime, then use an isothermal approximation
        Am(:) = (0.5_dp*(be(2:) + be(1:nlay)) * e0i(:))/be(1:nlay)
        Bm(:) = 0.0_dp
        Gp(:) = 0.0_dp
        Bp(:) = Am(:)
      elsewhere
        ! Use linear interpolants
        e1i(:) = del(:) - e0i(:)
        e1i_del(:) = e1i(:)/del(:) ! The equivalent to the linear in tau term

        Am(:) = e0i(:) - e1i_del(:) ! Am(k) = Gp(k), just indexed differently
        Bm(:) = e1i_del(:) ! Bm(k) = Bp(k), just indexed differently
        Gp(:) = Am(:)
        Bp(:) = Bm(:)
      end where

      !! Begin two-stream loops
      !! Perform downward loop first
      ! Top boundary condition - 0 flux downward from top boundary
      lw_down_g(1) = 0.0_dp
      do k = 1, nlay
        lw_down_g(k+1) = lw_down_g(k)*edel(k) + Am(k)*be(k) + Bm(k)*be(k+1) ! TS intensity
      end do

      !! Perform upward loop
      ! Lower boundary condition - internal heat definition Fint = F_down - F_up
      ! here we use the same condition but use intensity units to be consistent
      lw_up_g(nlev) = lw_down_g(nlev) + be_int
      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*edel(k) + Bp(k)*be(k) + Gp(k)*be(k+1) ! TS intensity
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      flx_down(:) = flx_down(:) + lw_down_g(:) * w(m)
      flx_up(:) = flx_up(:) + lw_up_g(:) * w(m)

    end do

    !! The flux is the intensity * pi
    flx_down(:) = pi * flx_down(:)
    flx_up(:) = pi * flx_up(:)

  end subroutine lw_shortchar_linear

  subroutine BB_integrate(n_b, Te, wn_e, be)
    implicit none

    integer, intent(in) :: n_b
    real(dp), intent(in) :: Te
    real(dp), dimension(n_b+1), intent(in) :: wn_e

    real(dp), dimension(n_b), intent(out) :: be

    integer :: ww, j, intitera
    real(dp), dimension(n_b+1) :: iB
    real(dp) :: x, x2, x3, itera, summ, dn

    !! Code for integrating the blckbody function between two wavenumbers
    !! This is a method that uses a sum convergence
    !! Taken from: spectralcalc.com/blackbody/inband_radiance.html

      if (Te < 1e-6_dp) then
        be(:) = 0.0_dp
        return
      end if

      do ww = 1, n_b+1

        x = c1 * 100.0_dp * wn_e(ww)/ Te
        x2 = x**2
        x3 = x**3

        itera = 2.0_dp + 20.0_dp/x
        if (itera > 150) then
          itera = 150
        end if
        intitera = int(itera)

        summ = 0.0_dp
        do j = 1, intitera + 1
          dn = 1.0_dp/real(j,kind=dp)
          summ = summ +  exp(-min(real(j,kind=dp)*x,300.0_dp))* &
          & (x3 + (3.0_dp * x2 + 6.0_dp*(x+dn)*dn)*dn)*dn
        end do

        iB(ww) = n2 * (Te/c1)**(4) * summ
      end do

      do ww = 1, n_b
        be(ww) = max(iB(ww+1) - iB(ww),0.0_dp)
      end do

  end subroutine BB_integrate

end module lw_sc_linear_mod