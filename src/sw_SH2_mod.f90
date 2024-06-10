!!!
! Elspeth KH Lee - Jun 2024: Initial version
! sw: Two stream spherical harmonic method with multiple scattering (Rooney et al. 2023) 
! following the PICASO method
!     Pros: 
!     Cons: 
!!!

!!! In development!!

module sw_SH2_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: fourpi = 4.0_dp * pi

  !! Use Two-Term HG function for sw
  logical, parameter :: TTHG = .False.

  private :: sw_SH_two_stream
  public :: sw_SH2

contains

  subroutine sw_SH2(nlay, nlev, nb, ng, gw, tau_e, mu_z, Finc, ssa, gg, a_surf, sw_up, sw_down, sw_net, asr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev, nb, ng
    real(dp), dimension(ng) :: gw
    real(dp), dimension(ng,nb,nlev), intent(in) :: tau_e
    real(dp), dimension(ng,nb,nlay), intent(in) :: ssa, gg
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
          call sw_SH_two_stream(nlay, nlev, Finc(b), mu_z(:), tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), a_surf(b), sw_down_g(g,:), sw_up_g(g,:))
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

  end subroutine sw_SH2

  subroutine sw_SH_two_stream(nlay, nlev, F0_in, mu_in, tau_in, w_in, g_in, w_surf_in, flx_down, flx_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0_in, w_surf_in
    real(dp), dimension(nlev), intent(in) :: tau_in, mu_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    !! Work variables
    integer :: l, k, i
    real(dp), dimension(nlev) :: cum_trans
    real(dp) :: mu_z
    integer, parameter :: nstr = 2
    real(dp), dimension(nstr) :: Pu0
    real(dp), dimension(nstr, nlay) :: a, bb, w_multi
    real(dp), dimension(nlay) :: dtau
    real(dp) :: surf_reflect, b_surface, b_top

    real(dp), dimension(2,nlay) :: eta
    real(dp), dimension(nlay) :: del, lam, expo, exptrm, q, Q1, Q2
    real(dp), dimension(nlay) :: Q1mn, Q2mn, Q1pl, Q2pl
    real(dp), dimension(nlay) :: zmn, zpl, zmn_up, zpl_up, zmn_down, zpl_down
    real(dp), dimension(nlev) :: expon, tau_e, dir

    real(dp), dimension(5,2*nlay) :: Mb
    real(dp), dimension(7,2*nlay) :: Mb_F
    real(dp), dimension(2*nlay) :: B, X

    real(dp), dimension(2*nlev, 2*nlay) :: F
    real(dp), dimension(2*nlev) :: G
    real(dp), dimension(2*nlay) :: F_bot
    real(dp) :: G_bot

    real(dp), dimension(2*nlev) :: flux_temp
    real(dp) :: flux_bot

    real(dp), parameter :: eps_20 = 1.0e-20_dp

    real(dp), dimension(nlay) :: sigma_sq, pmom2, fc, c
    real(dp), dimension(nlay) :: w0, hg
    integer :: info
    integer, dimension(2*nlay) :: ipiv

    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w_in(:) <= 1.0e-12_dp)) then

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

      return

    end if

    !! Calculate dtau in each layer
    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)

    !! Delta-M+ scaling (Following DISORT: Lin et al. 2018)
    !! Assume HG phase function for scaling (g**nstream)
    ! where (g_in(:) >= 1e-6_dp)
    !   fc(:) = g_in(:)**(nstr)
    !   pmom2(:) = g_in(:)**(nstr+1)
    !   sigma_sq(:) = real((nstr+1)**2 - nstr**2,dp) / &
    !   & ( log(fc(:)**2/pmom2(:)**2) )
    !   c(:) = exp(real(nstr**2,dp)/(2.0_dp*sigma_sq(:)))
    !   fc(:) = c(:)*fc(:)

    !   w0(:) = w_in(:)*((1.0_dp - fc(:))/(1.0_dp - fc(:)*w_in(:)))
    !   dtau(:) = (1.0_dp - w_in(:)*fc(:))*dtau(:)

    ! elsewhere
    !   w0(:) = w_in(:)
    !   fc(:) = 0.0_dp
    ! end where

    fc(:) = 0.0_dp

    hg(:) = g_in(:)
    w0(:) = w_in(:) 

    mu_z = mu_in(nlev) ! Can't do spherical geometry yet

    !! Reform edge optical depths
    tau_e(1) = 0.0_dp
    do k = 1, nlay
      tau_e(k+1) = tau_e(k) + dtau(k)
    end do

    dir(:) = F0_in * mu_z * exp(-tau_e(:)/mu_z)

    Pu0(1) = 1.0_dp
    Pu0(2) = -mu_z

    w_multi(1,:) = 1.0_dp
    w_multi(2,:) = 3.0_dp * (hg(:) - fc(:)) / (1.0_dp - fc(:))

    do l = 1, nstr
      a(l,:) = real(2*(l-1) + 1,dp) -  w0(:) * w_multi(l,:) + eps_20
      bb(l,:) = ((w0(:) * w_multi(l,:)) * 1.0_dp * Pu0(l)) / (4.0_dp*pi)
    end do

    surf_reflect = 0.0_dp

    b_surface = 0.0_dp + surf_reflect*mu_z*1.0_dp*exp(-tau_e(nlev)/mu_z)

    b_top = 0.0_dp

    lam(:) = sqrt(a(1,:)*a(2,:))
    expo(:) = min(lam(:)*dtau(:),35.0_dp)
    exptrm(:) = exp(-expo(:))

    del(:) = ((1.0_dp / mu_z)**2 - a(1,:)*a(2,:))
    eta(1,:) = (bb(2,:)/mu_z - a(2,:)*bb(1,:)) / del(:)
    eta(2,:) = (bb(1,:)/mu_z - a(1,:)*bb(2,:)) / del(:)

    q(:) = lam(:)/a(2,:)
    Q1(:) = (1.0_dp + 2.0_dp*q(:))*pi
    Q2(:) = (1.0_dp - 2.0_dp*q(:))*pi

    Q1mn(:) = Q1(:)*exptrm(:);  Q2mn(:) = Q2(:)*exptrm(:)
    Q1pl(:) = Q1(:)/exptrm(:);  Q2pl(:) = Q2(:)/exptrm(:)

    zmn(:) = (eta(1,:) - 2.0_dp*eta(2,:))*pi
    zpl(:) = (eta(1,:) + 2.0_dp*eta(2,:))*pi
    expon(:) = exp(-tau_e(:)/mu_z)
    zmn_up(:) = zmn(:) * expon(2:nlev)
    zpl_up(:) = zpl(:) * expon(2:nlev) 
    zmn_down(:) = zmn(:) * expon(1:nlay)
    zpl_down(:) = zpl(:) * expon(1:nlay) 

    Mb(:,:) = 0.0_dp
    B(:) = 0.0_dp

    Mb(3,1) = Q1(1)
    Mb(2,2) = Q2(1)
    B(1) = b_top - zmn_down(1)

    Mb(4, 2*nlay-1) = Q2mn(nlay) - surf_reflect*Q1mn(nlay)
    Mb(3, 2*nlay) = Q1pl(nlay) - surf_reflect*Q2pl(nlay)
    B(2*nlay) = b_surface - zpl_up(nlay) + surf_reflect*zmn_up(nlay)

    do i = 2, nlay
      Mb(1, 2*i) = -Q2(i)
      Mb(2, 2*i - 1) = -Q1(i)
      Mb(2, 2*i) = -Q1(i)
      Mb(3, 2*i - 1) = -Q2(i)
    end do

    do i = 1, nlay-1
      Mb(3, 2*i) = Q2pl(i)
      Mb(4, 2*i-1) = Q1mn(i)
      Mb(4, 2*i) = Q1pl(i)
      Mb(5, 2*i-1) = Q2mn(i)
    end do

    do i = 1, nlay-1
      B(2*i) = zmn_down(i+1) - zmn_up(i)
      B(2*i + 1) = zpl_down(i+1) - zpl_up(i)
    end do

    Mb_F(:,:) = 0.0_dp
    Mb_F(3,:) = Mb(1,:)
    Mb_F(4,:) = Mb(2,:)
    Mb_F(5,:) = Mb(3,:)  
    Mb_F(6,:) = Mb(4,:)
    Mb_F(7,:) = Mb(5,:)   

    call dgbsv(2*nlay, 2, 2, 1, Mb_F, 7, ipiv, B, 2*nlay, info)

    if (info == 0) then
        ! Success, B now contains the solution
    else
        print *, "An error occurred: ", info
    endif

    ! flux at bottom of atmosphere
    F_bot(2*nlay - 1) = Q2mn(nlay)
    F_bot(2*nlay) = Q1pl(nlay)
    G_bot = zpl_up(nlay)

    F(1, 1) = Q1(1)
    F(1, 2) = Q2(1)
    F(2, 1) = Q2(1)
    F(2, 2) = Q1(1)

    k = 0
    do i = 1, 2*nlay, 2
      F(i + 2, i) = Q1mn(k + 1)
      F(i + 2, i + 1) = Q2pl(k + 1)
      F(i + 3, i) = Q2mn(k + 1)
      F(i + 3, i + 1) = Q1pl(k + 1)
      k = k + 1
    end do

    G(1) = zmn_down(1)
    G(2) = zpl_down(1)

    do i = 1, nlay
      G(i*2 + 1) = zmn_up(i)
    end do

    do i = 2, nlay
      G(i*2) = zpl_up(i)
    end do

    flux_temp(:) = matmul(F(:,:),B(:)) + G(:)
    flux_bot = sum(F_bot(:)*X(:)) + G_bot

    do i = 1, nlev
      flx_down(i) = flux_temp(i*2-1)*F0_in 
      flx_up(i) = flux_temp(i*2)*F0_in
    end do

    ! Add direct beam
    flx_down(:) = flx_down(:) + dir(:)

  end subroutine sw_SH_two_stream

end module sw_SH2_mod