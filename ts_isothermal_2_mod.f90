!!!
! Elspeth KH Lee - May 2021
! Two-stream method following the isothermal layer approximation
! Pros: Very fast, better at high optical depths than 1st isothermal method
! Cons: Worse at low optical depth, worse at lower internal temperatures
!!!

module ts_isothermal_2_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: sb = 5.670374419e-8_dp

  public :: ts_isothermal_2
  private :: lw_grey_updown, sw_grey_down, linear_log_interp

contains

  subroutine ts_isothermal_2(nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, Beta_V, Beta_IR, net_F)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(3,nlev), intent(in) :: tau_Ve
    real(dp), dimension(2,nlev), intent(in) :: tau_IRe
    real(dp), dimension(3), intent(in) :: Beta_V
    real(dp), dimension(2), intent(in) :: Beta_IR
    real(dp), intent(in) :: F0, mu_z, Tint, AB

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: net_F

    !! Work variables
    integer :: i, b
    real(dp) :: Finc, be_int, Finc_b, be_int_b
    real(dp), dimension(nlev) :: Te, be, be_b
    real(dp), dimension(nlev) :: sw_down, sw_up, lw_down, lw_up
    real(dp), dimension(3,nlev) :: sw_down_b, sw_up_b
    real(dp), dimension(2,nlev) :: lw_down_b, lw_up_b
    real(dp), dimension(nlev) :: lw_net, sw_net

    !! Find temperature at layer edges through linear interpolation and extrapolation
    do i = 2, nlay
      call linear_log_interp(pe(i), pl(i-1), pl(i), Tl(i-1), Tl(i), Te(i))
      !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
    end do
    Te(1) = Tl(1) + (pe(1) - pe(2))/(pl(1) - pe(2)) * (Tl(1) - Te(2))
    Te(nlev) = Tl(nlay) + (pe(nlev) - pe(nlay))/(pl(nlay) - pe(nlay)) * (Tl(nlay) - Te(nlay))

    !! Shortwave flux calculation
    if (mu_z > 0.0_dp) then
      Finc = (1.0_dp - AB) * F0
      sw_down(:) = 0.0_dp
      do b = 1, 3
        Finc_b = Finc * Beta_V(b)
        call sw_grey_down(nlev, Finc_b, tau_Ve(b,:), mu_z, sw_down_b(b,:))
        sw_down(:) = sw_down(:) + sw_down_b(b,:)
      end do
    else
      sw_down(:) = 0.0_dp
    end if
    sw_up(:) = 0.0_dp ! sw_up is zero since we don't have shortwave scattering in this mode

    !! Longwave two-stream flux calculation
    be(:) = sb * Te(:)**4  ! Integrated planck function flux at levels
    be_int = sb * Tint**4 ! Integrated planck function flux for internal temperature
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp
    do b = 1, 2
      be_b(:) = be(:) * Beta_IR(b)
      be_int_b = be_int * Beta_IR(b)
      call lw_grey_updown(nlay, nlev, be_b, be_int_b, tau_IRe(b,:), lw_up_b(b,:), lw_down_b(b,:))
      lw_up(:) = lw_up(:) + lw_up_b(b,:)
      lw_down(:) = lw_down(:) + lw_down_b(b,:)
    end do

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

  end subroutine ts_isothermal_2

  subroutine lw_grey_updown(nlay, nlev, be, be_int, tau_IRe, lw_up, lw_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_IRe
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k
    real(dp), dimension(nlay) :: dtau, B0

    !! Prepare loop
    do k = 1, nlay
      dtau(k) = tau_IRe(k+1) - tau_IRe(k)
      B0(k) = 0.5_dp * (be(k+1) + be(k)) ! Layer average planck flux (isothermal approximation)
    end do

    !! First do the downward loop
    lw_down(1) = 0.0_dp
    do k = 1, nlay
       lw_down(k+1) = (2.0_dp * B0(k) * dtau(k))/(2.0_dp + dtau(k)) + &
       & lw_down(k) * (2.0_dp - dtau(k))/(2.0_dp + dtau(k))
    end do

    !! Perform upward loop
    ! Lower boundary condition - internal heat definition Fint = F_up - F_down
    lw_up(nlev) = lw_down(nlev) + be_int
    do k = nlay, 1, -1
      lw_up(k) = (2.0_dp * B0(k) * dtau(k))/(2.0_dp + dtau(k)) + &
      & lw_up(k+1) * (2.0_dp - dtau(k))/(2.0_dp + dtau(k))
    end do


  end subroutine lw_grey_updown

  subroutine sw_grey_down(nlev, Finc, tau_V, mu_z, sw_down)
    implicit none

    !! Input
    integer, intent(in) :: nlev
    real(dp), intent(in) :: Finc, mu_z
    real(dp), dimension(nlev), intent(in) :: tau_V

    !! Output
    real(dp), dimension(nlev), intent(out) :: sw_down

    sw_down(:) = Finc * mu_z * exp(-tau_V(:)/mu_z)

  end subroutine sw_grey_down

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    lxval = log10(xval)
    lx1 = log10(x1); lx2 = log10(x2)
    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / (lx2 - lx1)

    yval = 10.0_dp**((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)

  end subroutine linear_log_interp

end module ts_isothermal_2_mod
