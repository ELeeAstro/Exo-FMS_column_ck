!!!
! Elspeth KH Lee - Jun 2024 : Initial version
!                - Jul 2025 : n-stream DISORT 4.0.99, gas-giant lower boundary
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!
! lw: n-stream DISORT method - intended as an accuracy benchmark
!     Pros: Highly accurate reference solution, arbitrary stream number
!     Cons: Very slow compared to the two-stream methods
!
! Lower boundary
! --------------
! DISORT is written for a planetary surface, so its lower boundary emits
! sigma*BTEMP**4 and reflects with albedo ALBEDO. A gas giant has no surface -
! the requirement is instead that the interior supplies a prescribed NET
! upward flux,
!
!     F_up(nlev) - F_down(nlev) = sigma*Tint**4 (in-band: pi*B_band(Tint) for each band and g-ordinance),
!
! with downwelling radiation absorbed rather than reflected.
!
! We impose this on the *unmodified* DISORT by exploiting the linearity of the
! radiative transfer equation in its source terms. With ALBEDO = 0 the only
! boundary source is the lower-boundary emission x = B(BTEMP), and every flux
! is an affine function of it,
!
!     F(x) = P + Q*x.
!
! Two DISORT calls with different BTEMP therefore determine the whole family,
! and the combination
!
!     F = F1 + alpha*(F2 - F1),    alpha = (Fint - net1) / (net2 - net1)
!
! is the exact solution carrying the same atmospheric Planck sources (their
! weight is (1-alpha) + alpha = 1) and whatever lower-boundary emission is
! needed to hit the required net flux. No modification of DISORT is required,
! which matters for a routine whose purpose is to be a trustworthy reference.
!!!

module lw_disort_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4
  use BB_integrate_mod, only : BB_integrate
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

  !! Number of computational streams used by DISORT.
  !! Increase for a more accurate benchmark, at proportionally more cost.
  !! You will need to make clean and recompile if this is changed.
  integer, parameter :: nstr = 32

  !! Use the delta-M+ truncation (Lin et al. 2018). DISORT applies its own
  !! guards internally and falls back to plain delta-M where the moment
  !! sequence makes the delta-M+ fit unreliable (DISORT.f:482-497).
  logical, parameter :: use_deltamplus = .True.

  private
  public :: lw_disort

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
    integer :: b, g
    real(dp), dimension(nlev) :: Te
    real(dp), dimension(nlev) :: up1, dn1, up2, dn2
    real(dp), dimension(nb) :: be_int
    real(dp) :: Tbot1, Tbot2, net1, net2, alp, denom, Fint

    !! Use WENO4 method to (smoothly) interpolate layers to levels, in ln p
    Te(:) = interpolate_weno4(log(pe), log(pl), Tl, .False.)

    !! Edges are linearly extrapolated in log(p)-log(T) space.
    Te(1) = 10.0_dp**(log10(Te(2)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * &
      & log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Te(nlay)) + (log10(pe(nlev)/pe(nlay))/ &
      & log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    !! Integrated planck function intensity for internal temperature in each band
    call BB_integrate(nb, Tint, wn_e, be_int)

    !! Two lower-boundary emission temperatures spanning the solution family.
    !! They are scaled to the deepest atmospheric temperature so the pair stays
    !! well separated and well conditioned for any Tint, including Tint = 0.
    Tbot1 = Te(nlev)
    Tbot2 = Te(nlev) * 2.0_dp**0.25_dp

    !! Sum over each band and g-ordinance
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp
    do b = 1, nb

      !! Required in-band net internal flux
      Fint = pi * be_int(b)

      do g = 1, ng

        call call_DISORT(nlay, nlev, Te, tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), wn_e(b+1), wn_e(b), &
          & Tbot1, up1, dn1)
        call call_DISORT(nlay, nlev, Te, tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), wn_e(b+1), wn_e(b), &
          & Tbot2, up2, dn2)

        !! Solve for the lower-boundary emission that delivers the internal flux
        net1 = up1(nlev) - dn1(nlev)
        net2 = up2(nlev) - dn2(nlev)
        denom = net2 - net1

        !! A band with negligible emission at the lower boundary cannot be adjusted, use the first solve
        if (abs(denom) <= tiny(denom)) then
          alp = 0.0_dp
        else
          alp = (Fint - net1)/denom
        end if

        lw_up(:) = lw_up(:) + (up1(:) + alp*(up2(:) - up1(:))) * gw(g)
        lw_down(:) = lw_down(:) + (dn1(:) + alp*(dn2(:) - dn1(:))) * gw(g)

      end do
    end do

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_disort

  subroutine call_DISORT(nlay, nlev, Te, tau_e, ssa, gg, wn_lo, wn_hi, Tbot, flx_up, flx_down)
    implicit none

    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: Te, tau_e
    real(dp), dimension(nlay), intent(in) :: ssa, gg
    real(dp), intent(in) :: wn_lo, wn_hi
    real(dp), intent(in) :: Tbot
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! DISORT dimensioning
    integer, parameter :: nphi = 1
    integer, parameter :: nmom = nstr + 1
    integer :: numu, maxulv, i, k

    !! DISORT arguments - the DISORT toolchain is compiled with
    !! -fdefault-real-8, so its REAL declarations are 8 byte (see Makefile)
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
      !! Henyey-Greenstein phase function moments: pmom(k) = g**k
      pmom(0,i) = 1.0_dp
      do k = 1, nmom
        pmom(k,i) = real(gg(i)**k, dp)
      end do
    end do

    !! DISORT indexes temperature from 0 at the top of the column
    do i = 1, nlev
      temper(i-1) = real(Te(i), dp)
      utau(i) = real(tau_e(i), dp)
    end do

    !! Thermal source only - no solar beam in the longwave
    plank = .True.
    fbeam = 0.0_dp
    umu0 = 1.0_dp
    phi0 = 0.0_dp
    fisot = 0.0_dp
    phi(:) = 0.0_dp

    !! Band edges [cm-1] for the band-integrated Planck function
    wvnmlo = wn_lo
    wvnmhi = wn_hi

    !! Non-reflecting lower boundary emitting at Tbot. The physical gas-giant
    !! condition is imposed by the caller through the linear combination.
    lamber = .True.
    albedo = 0.0_dp
    btemp = real(Tbot, dp)

    !! No emission or reflection from the top boundary
    ttemp = 0.0_dp
    temis = 0.0_dp

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

    !! DISORT overwrites this argument internally (DISORT.f:474-497), so it
    !! must be reset on every call.
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

end module lw_disort_mod
