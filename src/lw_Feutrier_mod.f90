!!!
! Elspeth KH Lee - Jun 2024 : Initial version
! lw: Feutrier
!     Pros: 
!     Cons: 
!!!

module lw_AA_L_mod
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
  

  private :: lw_Feutrier_method, BB_integrate
  public :: lw_Feutrier

contains


  subroutine lw_Feutrier(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
        call lw_Feutrier_method(nlay, nlev, be(b,:), be_int(b), tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), &
          & lw_up_g(g,:), lw_down_g(g,:))
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

  end subroutine lw_Feutrier

  subroutine lw_Feutrier_method(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables - we follow the Hubeny indexing conventions
    integer :: d, l, k 
    real(dp), dimension(nlay) :: dtau, w0, hg
    integer, parameter :: n = 3

    real(dp), dimension(nlev,n,n) :: A, B, C
    real(dp), dimension(nlev,n) :: R, j
    real(dp), dimension(nlev) :: D, E

    w0(:) = w_in(:)
    hg(:) = g_in(:)

    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)

    A(:,:,:) = 0.0_dp
    B(:,:,:) = 0.0_dp
    C(:,:,:) = 0.0_dp
    R(:) = 0.0_dp

    do d = 1, nlev


      if (d == 1) then

        ! Special boundary conditions
        do l = 1, n
          R(k,l) = (1.0_dp - w0(k))*be(k)
          do k = 1, n

            if (k == l) then
              kdel = 1.0_dp
            else
              kdel = 0.0_dp 
            end if


            A(d,l,k) = 0.0_dp
            C(d,l,k) = 2.0_dp * (mu(l)/dtau_32)**2 * kdel
            B(d,l,k) = (1.0_dp + (2.0_dp*mu(l)/dtau_32) + 2.0_dp*(mu(l)/dtau_32)**2) * kdel - alph(d,l) * phi(d,k)

          end do
        end do

      else if (d == nlev) then

        ! Special boundary conditions
        do l = 1, n
          R(k,l) = (1.0_dp - w0(k))*be(k) + (2.0_dp*mu(l)/dtau_mh) * be_int
          do k = 1, n

            if (k == l) then
              kdel = 1.0_dp
            else
              kdel = 0.0_dp 
            end if


            A(d,l,k) = 2.0_dp * (mu(l)/dtau_mh)**2 * kdel
            C(d,l,k) = 0.0_dp
            B(d,l,k) = (1.0_dp + (2.0_dp*mu(l)/dtau_mh) + 2.0_dp*(mu(l)/dtau_mh)**2) * kdel - alph(d,l) * phi(d,k)

          end do
        end do

      else
        !! Find the values of the large tridiagonal matrix
        do l = 1, n
          R(k,l) = (1.0_dp - w0(k))*be(k)

          do k = 1, n

            if (k == l) then
              kdel = 1.0_dp
            else
              kdel = 0.0_dp 
            end if

            A(d,l,k) = mu(l) / (dtau(d-1)*dtau(d)) * kdel
            C(d,l,k) = mu(l) / (dtau(d+1)*dtau(d)) * kdel
            B(d,l,k) = kdel +  A(d,l,k) + C(d,l,k) - alph(d,l) * w(d,k)

          end do
        end do 
      end if

    end do


    do l = 1, n
      ! Find D and E for this mu and solve for j
      do d = 1, nlev

      end do
    end do



  end subroutine lw_Feutrier_method

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

end module lw_AA_L_mod