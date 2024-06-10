!!!
! Elspeth KH Lee - Jun 2024 : Initial version
!
! sw: Direct beam component only
!     Pros: Very very fast
!     Cons: Zero scattering
!!!

module sw_direct_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  private :: sw_direct_beam
  public :: sw_direct

contains

  subroutine sw_direct(nlev, nb, ng, gw, tau_e, mu_z, Finc, a_surf, sw_up, sw_down, sw_net, asr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlev, nb, ng
    real(dp), dimension(ng) :: gw
    real(dp), dimension(ng,nb,nlev), intent(in) :: tau_e
    real(dp), dimension(nlev) :: mu_z
    real(dp), dimension(nb) :: Finc, a_surf

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down, sw_net
    real(dp), intent(out) :: asr

    !! Work variables
    integer :: b, g
    real(dp), dimension(nb,nlev) :: sw_down_b, sw_up_b
    real(dp), dimension(ng,nlev) :: sw_down_g, sw_up_g


    !! Shortwave flux calculation
    if (mu_z(nlev) > 0.0_dp) then
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
      do b = 1, nb
        sw_down_b(b,:) = 0.0_dp
        sw_up_b(b,:) = 0.0_dp
        if (Finc(b) <= 0.0_dp) then
          cycle
        end if
        do g = 1, ng
          call sw_direct_beam(nlev, Finc(b), mu_z(:), tau_e(g,b,:), a_surf(b), sw_down_g(g,:), sw_up_g(g,:))
          sw_down_b(b,:) = sw_down_b(b,:) + sw_down_g(g,:) * gw(g)
          sw_up_b(b,:) = sw_up_b(b,:) + sw_up_g(g,:) * gw(g)
        end do
        sw_down(:) = sw_down(:) + sw_down_b(b,:)
        sw_up(:) = sw_up(:) + sw_up_b(b,:)
      end do
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Net sw flux
    sw_net(:) = sw_up(:) - sw_down(:)

    !! Absorbed Stellar Radiation (ASR)
    asr = sw_down(1) - sw_up(1)

  end subroutine sw_direct

  subroutine sw_direct_beam(nlev, F0_in, mu_in, tau_in, w_surf_in, flx_down, flx_up)
    implicit none

    integer, intent(in) :: nlev
    real(dp), intent(in) :: F0_in, w_surf_in
    real(dp), dimension(nlev), intent(in) :: tau_in, mu_in

    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    integer :: k
    real(dp), dimension(nlev) :: cum_trans

    !! Calculate direct beam component only

    ! Check mu_z for spherical corrections
    if (mu_in(nlev) == mu_in(1)) then
      ! No zenith correction, use regular method
      flx_down(:) = F0_in* mu_in(nlev) * exp(-tau_in(:)/mu_in(nlev))
    else
      ! Zenith angle correction, use cumulative transmission function
      cum_trans(1) = tau_in(1)/mu_in(1)
      do k = 1, nlev-1
        cum_trans(k+1) = cum_trans(k) + (tau_in(k+1) - tau_in(k))/mu_in(k+1)
      end do
      do k = 1, nlev
        flx_down(k) = F0_in * mu_in(nlev) * exp(-cum_trans(k))
      end do
    end if

    flx_down(nlev) = flx_down(nlev) * (1.0_dp - w_surf_in) ! The surface flux for surface heating is the amount of flux absorbed by surface
    flx_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo


  end subroutine sw_direct_beam


end module sw_direct_mod