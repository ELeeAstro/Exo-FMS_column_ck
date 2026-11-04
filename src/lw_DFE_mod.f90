!!!
! Elspeth KH Lee - July 2024 : Initial version
! lw: Absorption-only discontinuous finite-element method.
!     Scattering is not represented by this solver.
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!                - Nov 2026 : Optimisations - reciprocals and source terms outside the sweep recurrences,
!                             contiguous per g-ordinance inputs
!!!

module lw_DFE_mod
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
  real(dp), parameter :: dtau_floor = 1.0e-12_dp

  !! Gauss quadrature variables, cosine angle values (uarr) and weights (w)
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these
  
  !! Optimised quadrature for 1 node - 2 stream (Hogan 2024)
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.6096748751_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)

  !! Gauss–Jacobi-5 quadrature for 2 nodes - 4 stream (Hogan 2024)
  !integer, parameter :: nmu = 2
  !real(dp), dimension(nmu), parameter :: uarr = (/0.2509907356_dp, 0.7908473981_dp/)
  !real(dp), dimension(nmu), parameter :: w = (/0.2300253764_dp, 0.7699746236_dp/)

  !! Gauss–Jacobi-5 quadrature for 3 nodes - 6 stream (Hogan 2024)
  ! integer, parameter :: nmu = 3
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.1024922169_dp, 0.4417960320_dp, 0.8633751621_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/0.0437820218_dp, 0.3875796738_dp, 0.5686383044_dp/)

  !! Gauss–Jacobi-5 quadrature for 4 nodes - 8 stream (Hogan 2024)
  integer, parameter :: nmu = 4
  real(dp), dimension(nmu), parameter :: uarr = (/0.0454586727_dp, 0.2322334416_dp, 0.5740198775_dp, 0.9030775973_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.0092068785_dp, 0.1285704278_dp, 0.4323381850_dp, 0.4298845087_dp /)
  real(dp), dimension(nmu), parameter :: ruarr = 1.0_dp/uarr

  private :: lw_discontinuous_finite_element
  public :: lw_DFE

contains


  subroutine lw_DFE(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, Tint, lw_up, lw_down, lw_net, olr)
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
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g
    real(dp), dimension(nlev) :: be_b, tau_g

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
        ! Contiguous copy (avoids an array temporary for the strided slice)
        tau_g(:) = tau_e(g,b,:)
        call lw_discontinuous_finite_element(nlay, nlev, be_b, be_int(b), tau_g, lw_up_g(:), lw_down_g(:))
        lw_up(:) = lw_up(:) + lw_up_g(:) * gw(g)
        lw_down(:) = lw_down(:) + lw_down_g(:) * gw(g)
      end do
    end do

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_DFE

  subroutine lw_discontinuous_finite_element(nlay, nlev, be, be_int, tau_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables and arrays
    integer :: k, m
    real(dp), dimension(nlay) :: dtau
    real(dp), dimension(nlay) :: del, ra, b, sdn, sup, pc
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g, Ip, Im
    real(dp), dimension(2:nlay) :: wk

    !! Calculate dtau in each layer
    dtau(:) = tau_in(2:) - tau_in(1:nlay)
    where (dtau(:) >= 0.0_dp .and. dtau(:) < dtau_floor)
      dtau(:) = dtau_floor
    end where

    ! Zero the total flux arrays
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nmu

      del(:) = dtau(:)*ruarr(m)

      !! Find 1/a and b, and the layer source terms, outside the recurrences
      b(:) = del(:)*(del(:) + 1.0_dp)
      ra(:) = 1.0_dp/(del(:)**2 + 2.0_dp*del(:) + 2.0_dp)
      pc(:) = 2.0_dp*(del(:) + 1.0_dp)*ra(:)

      !! Interface weights of the incoming-side intensity, del_k/(del_k + del_k-1)
      do k = 2, nlay
        if (abs(del(k) + del(k-1)) > tiny(1.0_dp)) then
          wk(k) = del(k)/(del(k) + del(k-1))
        else
          wk(k) = 1.0_dp
        end if
      end do

      !! Perform the Ip and Im calculation for the downward stream
      sdn(:) = (del(:)*be(1:nlay) + b(:)*be(2:nlev))*ra(:)
      sup(:) = (b(:)*be(1:nlay) - del(:)*be(2:nlev))*ra(:)
      Im(1) = 0.0_dp ! No external radiation
      do k = 1, nlay
        Im(k+1) = 2.0_dp*ra(k)*Im(k) + sdn(k)
      end do
      Ip(1:nlay) = pc(:)*Im(1:nlay) + sup(:)

      lw_down_g(1) = Im(1)
      do k = 2, nlay
        lw_down_g(k) = Im(k)*wk(k) + Ip(k)*(1.0_dp - wk(k))
      end do
      lw_down_g(nlev) = Im(nlev)

      !! Perform the Ip and Im calculation for the upward stream.
      !! Gas-giant lower boundary: the net upward intensity is B(Tint).
      sdn(:) = (del(:)*be(2:nlev) + b(:)*be(1:nlay))*ra(:)
      sup(:) = (b(:)*be(2:nlev) - del(:)*be(1:nlay))*ra(:)
      Im(nlev) = lw_down_g(nlev) + be_int
      do k = nlev, 2, -1
        Im(k-1) = 2.0_dp*ra(k-1)*Im(k) + sdn(k-1)
      end do
      Ip(2:nlev) = pc(:)*Im(2:nlev) + sup(:)

      lw_up_g(1) = Im(1)
      do k = 2, nlay
        lw_up_g(k) = Im(k)*wk(k) + Ip(k)*(1.0_dp - wk(k))
      end do
      lw_up_g(nlev) = Im(nlev)
      
      
      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      flx_down(:) = flx_down(:) + lw_down_g(:) * w(m)
      flx_up(:) = flx_up(:) + lw_up_g(:) * w(m)

    end do

    !! The flux is the intensity * pi
    flx_down(:) = pi * flx_down(:)
    flx_up(:) = pi * flx_up(:)

  end subroutine lw_discontinuous_finite_element

end module lw_DFE_mod
