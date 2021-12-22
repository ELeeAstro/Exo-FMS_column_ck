!!!
! Elspeth KH Lee - Sep 2020
! Module for dry convective adjustment schemes
! Currently has Ray Pierrehumbert's python notebook implementation
!
!!

module dry_conv_adj_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  integer, parameter :: itermax = 5
  real(dp), parameter :: small = 1.0e-6_dp

contains

  subroutine Ray_dry_adj(nlay, nlay1, t_step, kappa, Tl, pl, pe, dT_conv)
    implicit none

    integer, intent(in) :: nlay, nlay1
    real(dp), intent(in) :: kappa, t_step
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlay1), intent(in) :: pe
    real(dp), dimension(nlay), intent(out) :: dT_conv

    integer :: i, iter
    logical :: did_adj

    real(dp), dimension(nlay) :: Tl_cc, d_p
    real(dp) :: pfact, Tbar

    Tl_cc = Tl

    do i = 1, nlay
      d_p(i) = pe(i+1) - pe(i)
    end do

    do iter = 1, itermax

      did_adj = .False.

      !! Downward pass
      do i = 1, nlay-1
        pfact = (pl(i)/pl(i+1))**kappa
        if (Tl_cc(i) < (Tl_cc(i+1)*pfact - small)) then
          Tbar = (d_p(i)*Tl_cc(i) + d_p(i+1)*Tl_cc(i+1)) &
             & / (d_p(i) + d_p(i+1))
          Tl_cc(i+1) = (d_p(i) + d_p(i+1))*Tbar &
             & / (d_p(i+1)+pfact*d_p(i))
          Tl_cc(i) = Tl_cc(i+1) * pfact
          did_adj = .True.
        end if
      end do

      !! Upward pass
      do i = nlay-1, 1, -1
        pfact = (pl(i)/pl(i+1))**kappa
        if (Tl_cc(i) < (Tl_cc(i+1)*pfact - small)) then
          Tbar = (d_p(i)*Tl_cc(i) + d_p(i+1)*Tl_cc(i+1)) &
             & / (d_p(i) + d_p(i+1))
          Tl_cc(i+1) = (d_p(i) + d_p(i+1))*Tbar &
             & / (d_p(i+1)+pfact*d_p(i))
          Tl_cc(i) = Tl_cc(i+1) * pfact
          did_adj = .True.
        end if
      end do

      !! If no adjustment required, exit the loop
      if (did_adj .eqv. .False.) then
        exit
      end if

    end do

    !! Change in temperature is Tl_cc - Tl
    !! adjust on timescale of 1 timestep (i.e. instant correction across one timestep)
    dT_conv(:) = (Tl_cc(:) - Tl(:))/t_step

  end subroutine Ray_dry_adj

end module dry_conv_adj_mod
