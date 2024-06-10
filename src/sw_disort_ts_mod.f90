!!!
! Elspeth KH Lee - Jun 2024
! sw: two-stream disort method with modifications from Xianyu Tan.
!     Pros: Accurate and reasonably fast multiple scattering method
!     Cons: Two-stream only, can be inaccurate
!!!

module sw_disort_ts_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  private
  public :: sw_disort_ts

contains

  subroutine sw_disort_ts(nlay, nlev, nb, ng, gw, tau_e, mu_z, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
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
    integer :: b, g, i

    !! Conversion arrays from FMS to DISORT dimensions and work variables
    integer, parameter :: maxcly=200, maxulv=201
    real(dp), dimension(0:maxcly) :: Te_0
    real(dp) :: wvnmlo, wvnmhi
    real(dp), dimension(maxcly) :: dtauc, utau
    real(dp), dimension(maxcly) :: ggg, ssalb
    real(dp), dimension(maxulv) :: sw_net_d
    real(dp), dimension(nb,maxulv) :: sw_net_b
    real(dp), dimension(ng,maxulv) :: sw_net_g
    real(dp) :: umu0, fbeam, olr_g, Tint
    logical :: planck

    !! Shortwave flux calculation
    if (mu_z(nlev) > 0.0_dp) then
      planck = .False.
      umu0 =  mu_z(nlev)
      sw_net_d(:) = 0.0_dp
      !! Initalise arrays
      ggg(:) = 0.0_dp
      ssalb(:) = 0.0_dp
      utau(:) = 0.0_dp
      dtauc(:) = 0.0_dp
      Te_0(:) = 0.0_dp
      do b = 1, nb
        if (Finc(b) <= 0.0_dp) then
          cycle
        end if
        sw_net_b(b,:) = 0.0_dp
        fbeam = Finc(b)
        wvnmlo = 0.0_dp
        wvnmhi = 0.0_dp
        Tint = 0.0_dp
        do g = 1, ng
          ggg(1:nlay) = gg(g,b,:)
          ssalb(1:nlay) = ssa(g,b,:)
          utau(1:nlev) = tau_e(g,b,:)
          do i = 1, nlay
            dtauc(i) = (tau_e(g,b,i+1) - tau_e(g,b,i))
          end do
          call CALL_TWOSTR (nlay,Te_0,ggg,ssalb,dtauc,nlev,utau,planck,wvnmlo,wvnmhi,Tint,fbeam,umu0,sw_net_g(g,:),olr_g)
          sw_net_b(b,:) = sw_net_b(b,:) + sw_net_g(g,:) * gw(g)
        end do
        sw_net_d(:) = sw_net_d(:) + sw_net_b(b,:)
      end do
    else
      sw_net_d(:) = 0.0_dp
    end if

    sw_net(:) = sw_net_d(1:nlev)
    sw_up(:) = 0.0_dp
    sw_down(:) = 0.0_dp

    !! Net sw flux
    !sw_net(:) = sw_up(:) - sw_down(:)

    !! Absorbed Stellar Radiation (ASR)
    asr = 0.0_dp !sw_down(1) - sw_up(1)

  end subroutine sw_disort_ts

end module sw_disort_ts_mod