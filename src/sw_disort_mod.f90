!!!
! Elspeth KH Lee - Jun 2024 : Initial version
!                - Jul 2025 : n-stream DISORT 4.0.99
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!
! sw: n-stream DISORT method - intended as an accuracy benchmark
!     Pros: Highly accurate reference solution, arbitrary stream number
!     Cons: Very slow, single zenith angle so can't follow the zcorr beam path (unlike sw_Toon89 and sw_adding)
!
! Gas giant lower boundary - albedo = 0, so nothing is reflected back up (sw_up(nlev) = 0)
!!!

module sw_disort_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Number of DISORT streams - more is more accurate but slower
  !! make sure to make clean and recompile if you change this
  integer, parameter :: nstr = 32

  !! Use delta-M+ scaling (Lin et al. 2018) - DISORT falls back to delta-M itself where delta-M+ is unreliable
  logical, parameter :: use_deltamplus = .True.

  private
  public :: sw_disort

contains

  subroutine sw_disort(nlay, nlev, nb, ng, gw, tau_e, mu_z, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
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
    integer :: b, g
    real(dp), dimension(nlev) :: sw_up_g, sw_down_g

    sw_up(:) = 0.0_dp
    sw_down(:) = 0.0_dp

    !! Shortwave flux calculation - sum over each band and g-ordinance
    if (mu_z(nlev) > 0.0_dp) then
      do b = 1, nb
        if (Finc(b) <= 0.0_dp) cycle
        do g = 1, ng
          call call_DISORT(nlay, nlev, tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), mu_z(nlev), Finc(b), &
            & sw_up_g(:), sw_down_g(:))
          sw_up(:) = sw_up(:) + sw_up_g(:) * gw(g)
          sw_down(:) = sw_down(:) + sw_down_g(:) * gw(g)
        end do
      end do
    end if

    !! Net sw flux
    sw_net(:) = sw_up(:) - sw_down(:)

    !! Absorbed Stellar Radiation (ASR)
    asr = sw_down(1) - sw_up(1)

  end subroutine sw_disort

  subroutine call_DISORT(nlay, nlev, tau_e, ssa, gg, mu_in, F0_in, flx_up, flx_down)
    implicit none

    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: tau_e
    real(dp), dimension(nlay), intent(in) :: ssa, gg
    real(dp), intent(in) :: mu_in, F0_in
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! DISORT dimensioning
    integer, parameter :: nphi = 1
    integer, parameter :: nmom = nstr + 1
    integer :: numu, maxulv, i, k

    !! DISORT arguments - DISORT is compiled with -fdefault-real-8 (see Makefile)
    logical :: usrang, usrtau, onlyfl, plank, lamber, deltamplus, do_pseudo_sphere
    logical :: prnt(5)
    integer :: ibcnd
    real(dp) :: dtauc(nlay), ssalb(nlay), pmom(0:nmom,nlay), temper(0:nlay)
    real(dp) :: utau(nlev), h_lyr(0:nlay)
    real(dp) :: wvnmlo, wvnmhi, umu0, phi0, fbeam, fisot, albedo, btemp, ttemp, temis
    real(dp) :: earth_radius, accur
    real(dp) :: phi(nphi), umu(nstr)
    real(dp) :: rfldir(nlev), rfldn(nlev), flup(nlev), dfdt(nlev), uavg(nlev)
    real(dp) :: uu(nstr,nlev,nphi), albmed(nstr), trnmed(nstr)
    real(dp) :: rhoq(nstr/2,0:nstr/2,0:nstr-1), rhou(nstr,0:nstr/2,0:nstr-1)
    real(dp) :: emust(nstr), bemst(nstr/2), rho_accurate(nstr,nphi)
    character(len=127) :: header

    numu = nstr
    maxulv = nlev

    !! Layer optical properties
    do i = 1, nlay
      dtauc(i) = real(tau_e(i+1) - tau_e(i), dp)
      ssalb(i) = real(ssa(i), dp)
      ! Henyey-Greenstein phase function moments, pmom(k) = g**k
      pmom(0,i) = 1.0_dp
      do k = 1, nmom
        pmom(k,i) = real(gg(i)**k, dp)
      end do
    end do

    do i = 1, nlev
      utau(i) = real(tau_e(i), dp)
    end do

    !! Stellar beam only - no thermal source in the shortwave
    plank = .False.
    temper(:) = 0.0_dp
    wvnmlo = 0.0_dp
    wvnmhi = 0.0_dp
    btemp = 0.0_dp
    ttemp = 0.0_dp
    temis = 0.0_dp

    !! FBEAM is the flux normal to the beam (top of atmosphere flux is FBEAM*UMU0)
    fbeam = real(F0_in, dp)
    umu0 = real(mu_in, dp)
    phi0 = 0.0_dp
    fisot = 0.0_dp
    phi(:) = 0.0_dp

    !! Gas giant lower boundary - no reflection
    lamber = .True.
    albedo = 0.0_dp

    usrang = .False.
    usrtau = .False.
    onlyfl = .True.
    ibcnd = 0
    do_pseudo_sphere = .False.
    earth_radius = 6371.0_dp
    h_lyr(:) = 0.0_dp
    accur = 0.0_dp
    header = ''
    prnt(:) = .False.

    !! DISORT overwrites deltamplus internally, so reset it every call
    deltamplus = use_deltamplus

    rfldir(:) = 0.0_dp ; rfldn(:) = 0.0_dp ; flup(:) = 0.0_dp
    dfdt(:) = 0.0_dp ; uavg(:) = 0.0_dp ; uu(:,:,:) = 0.0_dp
    rhoq(:,:,:) = 0.0_dp ; rhou(:,:,:) = 0.0_dp
    emust(:) = 0.0_dp ; bemst(:) = 0.0_dp ; rho_accurate(:,:) = 0.0_dp
    albmed(:) = 0.0_dp ; trnmed(:) = 0.0_dp ; umu(:) = 0.0_dp

    call DISORT(nlay, nmom, nstr, numu, nphi, maxulv, &
      & usrang, usrtau, ibcnd, onlyfl, prnt, &
      & plank, lamber, deltamplus, do_pseudo_sphere, &
      & dtauc, ssalb, pmom, temper, wvnmlo, wvnmhi, &
      & utau, umu0, phi0, umu, phi, fbeam, &
      & fisot, albedo, btemp, ttemp, temis, &
      & earth_radius, h_lyr, &
      & rhoq, rhou, rho_accurate, bemst, emust, &
      & accur, header, &
      & rfldir, rfldn, flup, dfdt, uavg, uu, &
      & albmed, trnmed)

    do i = 1, nlev
      flx_up(i) = real(flup(i), dp)
      flx_down(i) = real(rfldir(i), dp) + real(rfldn(i), dp)
    end do

  end subroutine call_DISORT

end module sw_disort_mod
