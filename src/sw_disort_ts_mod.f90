!!!
! Elspeth KH Lee - Jun 2024 : Initial version
!                - Nov 2026 : Pass a local copy of nlev to TWOSTR
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!
! sw: Two-stream DISORT method with modifications from Xianyu Tan
!     Pros: Multiple scattering method, independent check of the sw_Toon89 scheme (same answer)
!     Cons: Two-stream only, ~3x slower than sw_Toon89
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
    real(dp), intent(in) :: mu_z
    real(dp), dimension(nb), intent(in) :: Finc

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down, sw_net
    real(dp), intent(out) :: asr

    !! Work variables
    integer :: b, g, i, ntau

    !! Conversion arrays from FMS to DISORT dimensions and work variables
    integer, parameter :: maxcly=200, maxulv=201
    real(dp), dimension(0:maxcly) :: Te_0
    real(dp) :: wvnmlo, wvnmhi
    real(dp), dimension(maxcly) :: dtauc
    real(dp), dimension(maxcly) :: ggg, ssalb
    real(dp), dimension(maxulv) :: sw_net_d, sw_up_d, sw_down_d, utau
    real(dp), dimension(maxulv) :: sw_net_g, sw_up_g, sw_down_g
    real(dp) :: umu0, fbeam, Tint
    logical :: planck

    sw_net_d(:) = 0.0_dp
    sw_up_d(:) = 0.0_dp
    sw_down_d(:) = 0.0_dp

    !! Shortwave flux calculation - sum over each band and g-ordinance
    if (mu_z > 0.0_dp) then

      planck = .False.
      umu0 =  mu_z

      !! Initialise arrays
      ggg(:) = 0.0_dp
      ssalb(:) = 0.0_dp
      utau(:) = 0.0_dp
      dtauc(:) = 0.0_dp
      Te_0(:) = 0.0_dp

      wvnmlo = 0.0_dp
      wvnmhi = 0.0_dp
      Tint = 0.0_dp

      do b = 1, nb
        if (Finc(b) <= 0.0_dp) cycle
        fbeam = Finc(b)
        do g = 1, ng
          ggg(1:nlay) = gg(g,b,:)
          ssalb(1:nlay) = ssa(g,b,:)
          utau(1:nlev) = tau_e(g,b,:) - tau_e(g,b,1)
          do i = 1, nlay
            dtauc(i) = (tau_e(g,b,i+1) - tau_e(g,b,i))
          end do
          ! Gas giant lower boundary (albedo = 0), TWOSTR overwrites ntau so we pass a copy of nlev
          ntau = nlev
          call call_twostr(nlay,0.0_dp,Te_0,ggg,ssalb,dtauc,ntau,utau,planck,wvnmlo,wvnmhi,Tint,fbeam,umu0, &
            & sw_net_g(:),sw_up_g(:),sw_down_g(:))
          sw_net_d(:) = sw_net_d(:) + sw_net_g(:) * gw(g)
          sw_up_d(:) = sw_up_d(:) + sw_up_g(:) * gw(g)
          sw_down_d(:) = sw_down_d(:) + sw_down_g(:) * gw(g)
        end do
      end do

    end if

    sw_net(:) = sw_net_d(1:nlev)
    sw_up(:) = sw_up_d(1:nlev)
    sw_down(:) = sw_down_d(1:nlev)

    !! Absorbed Stellar Radiation (ASR)
    asr = sw_down(1) - sw_up(1)

  end subroutine sw_disort_ts

end module sw_disort_ts_mod
