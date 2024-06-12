!!!
! Elspeth KH Lee - Jun 2024
! sw: n-stream disort method 
!     Pros: Highly Accurate, benchmark method
!     Cons: Extreamly slow - Can't do internal temperature properly (needs modification) - so just use net method
!!!

module lw_disort_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4  
  implicit none

  !! Precision variables
  integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64

  real(dp), parameter :: sb = 5.670374419e-8_dp

  private :: call_DISORT
  public :: lw_DISORT

contains


  subroutine lw_disort(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
    integer :: b, g, i
    real(dp), dimension(nlev) :: Te

    !! Conversion arrays from FMS to DISORT dimensions and work variables
    real(sp), dimension(0:nlay) :: Te_0
    real(sp) :: wvnmlo, wvnmhi
    real(sp), dimension(nlay) :: dtauc
    real(sp), dimension(nlev) :: utau
    real(sp), dimension(nlay) :: ggg, ssalb
    real(dp), dimension(nlev) :: lw_net_d
    real(dp), dimension(nb,nlev) :: lw_net_b
    real(dp), dimension(ng,nlev) :: lw_net_g
    real(sp) :: umu0, fbeam, Tint_sp
    real(dp) :: olr_g
    logical :: planck


    !! Use WENO4 method to (smoothly) interpolate layers to levels
    Te(:) = interpolate_weno4(pe, pl, Tl, .False.)

    !! Edges are linearly interpolated to avoid overshoot
    Te(1) = 10.0_dp**(log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Tl(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    Te_0(:) = 0.0_dp
    Te_0(0:nlay) = real(Te(1:nlev))

    !! Longwave two-stream flux calculation
    planck = .True.
    fbeam = 0.0_sp
    umu0 = 0.0_sp
    lw_net_d(:) = 0.0_dp
    olr = 0.0_dp
    ggg(:) = 0.0_sp
    ssalb(:) = 0.0_sp
    utau(:) = 0.0_sp
    dtauc(:) = 0.0_sp
    Tint_sp = Tint
    do b = 1, nb
      lw_net_b(b,:) = 0.0_dp
      wvnmlo = real(wn_e(b+1))
      wvnmhi = real(wn_e(b))
      do g = 1, ng
        do i = 1, nlay
          ggg(i) = real(gg(g,b,i))
          ssalb(i) = real(ssa(g,b,i))
          dtauc(i) = real((tau_e(g,b,i+1) - tau_e(g,b,i)))
        end do
        do i = 1, nlev
          utau(i) = real(tau_e(g,b,i))
        end do
        call call_DISORT(nlay,real(Te),ggg,ssalb,dtauc,nlev,utau,planck,wvnmlo,wvnmhi,real(Tint),fbeam,umu0,lw_net_g(g,:),olr_g)
        lw_net_b(b,:) = lw_net_b(b,:) + lw_net_g(g,:) * gw(g)
        olr = olr + olr_g * gw(g)
      end do
      lw_net_d(:) = lw_net_d(:) + lw_net_b(b,:)
    end do

    lw_net(:) = lw_net_d(1:nlev)
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp

    !! Apply net internal flux at lower boundary
    !lw_net(nlev) = sb * Tint**4

  end subroutine lw_disort

  subroutine call_DISORT(nlyr,temper,gg,ssalb,dtauc,ntau,utau, &
     &     plank,wvnmlo,wvnmhi,Tint,fbeam,umu0,fnetup,olr)
!     ------------------ INPUT VARIABLES -----------------
      LOGICAL plank
      INTEGER nlyr,ntau
      REAL*4 temper(0:nlyr),utau(ntau)
      REAL*4 gg(nlyr),ssalb(nlyr),dtauc(nlyr)
      REAL*4 fbeam,umu0,Tint,wvnmlo,wvnmhi
!     ------------------ OUTPUT VARIABLES ----------------
      REAL*8 fnetup(ntau), olr

!     ------------------ LOCAL VARIABLES -----------------
      LOGICAL :: deltamplus = .true.
      LOGICAL :: usrtau = .false.
      logical :: usrang = .false.
      logical :: lamber = .true.
      logical :: onlyfl = .true.
      logical :: do_pseudo_sphere = .false.
      LOGICAL prnt(5)

      REAL*4 :: albedo = 1.0_sp  ! 1 if use net upward flux
      REAL*4 :: fisot = 0.0_sp
      REAL*4 :: temis = 0.0_sp
      REAL*4 :: ttemp = 0.0_sp

      integer, parameter :: nphi = 1
      REAL*4 :: phi0 = 0.0_sp
      real*4 phi(nphi)

      integer :: ibcnd = 0

      integer, parameter :: nstr = 24
      integer, parameter :: nmom = nstr+1
      integer, parameter :: numu = nstr

      REAL*4 btemp,radius
      REAL*4 rfldir(ntau), rfldn(ntau), flup(ntau), &
     &     dfdt(ntau), uavg(ntau), uu(numu, ntau, nphi), &
     &     albmed(numu), trnmed(numu)
      REAL*4 taucum
      INTEGER i,j,k
      CHARACTER header*127

      real*4 pmom(0:nmom,nlyr), h_lyr(0:nlyr)
      real*4 rhoq(nstr/2, 0:nstr/2, 0:(nstr-1)), rhou(numu, 0:nstr/2, 0:(nstr-1))
      real*4 emust(numu), bemst(nstr/2), rho_accurate(numu, nphi)
      real*4 umu(numu)
      real*4 :: accur = 0.001_dp

      real*4 :: earth_radius = 6371.0_dp

      !! Find phase function Legendre expansion coefficients for each layer, assume HG function
      do i = 1, nlyr
        call GETMOM(3, gg(i), nmom, pmom(:,i) )
      end do

!     INITIALIZATION FOR FNETUP
      DO i = 1,nlyr+1
         fnetup(i) = 0.0_dp
      ENDDO

      taucum = 0.0_dp
      do i = 1,nlyr
         taucum = taucum + dtauc(i)
      enddo
      if (utau(ntau).gt.taucum) utau(ntau) = taucum
      btemp = Tint
      prnt(1) = .false.
      prnt(2) = .false.
      prnt(3) = .false.
      prnt(4) = .false.
      prnt(5) = .false.
      header = ''

    CALL DISORT( nlyr, NMOM, NSTR, NUMU, NPHI, ntau,           &
              & USRANG, USRTAU, IBCND, ONLYFL, PRNT,          &
              &   PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE,  &
              &   DTAUC, SSALB, PMOM, TEMPER, WVNMLO, WVNMHI,   &
              &   UTAU, UMU0, PHI0, UMU, PHI, FBEAM,            &
              &   FISOT, ALBEDO, BTEMP, TTEMP, TEMIS,           &
              &   EARTH_RADIUS, H_LYR,                          &
              &   RHOQ, RHOU, RHO_ACCURATE, BEMST, EMUST,       &
              &   ACCUR,  HEADER,                               &
              &   RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,          &
              &   ALBMED, TRNMED )

      do i = 1, ntau
         fnetup(i) = flup(i) - rfldir(i) - rfldn(i)
         !print*, i, flup(i), rfldir(i), rfldn(i)
      enddo

      olr = flup(1)

    end subroutine call_DISORT

end module lw_disort_mod