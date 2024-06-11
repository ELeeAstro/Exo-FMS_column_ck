!!!
! Elspeth KH Lee - Jun 2024: Initial version
! sw: Four stream spherical harmonic method with multiple scattering (Rooney et al. 2023) 
! following the PICASO method
!     Pros: 
!     Cons: 
!!!

!!! In development!!

module sw_SH4_mod
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

  private :: sw_SH_four_stream
  public :: sw_SH4

contains

  subroutine sw_SH4(nlay, nlev, nb, ng, gw, tau_e, mu_z, Finc, ssa, gg, a_surf, sw_up, sw_down, sw_net, asr)
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
          call sw_SH_four_stream(nlay, nlev, Finc(b), mu_z(:), tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), a_surf(b), & 
            & sw_down_g(g,:), sw_up_g(g,:))
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

  end subroutine sw_SH4

  subroutine sw_SH_four_stream(nlay, nlev, F0_in, mu_in, tau_in, w_in, g_in, w_surf_in, flx_down, flx_up)
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
    integer, parameter :: nstr = 4
    real(dp), dimension(nstr) :: Pu0
    real(dp), dimension(nstr, nlay) :: a, bb, w_multi
    real(dp), dimension(nlay) :: dtau
    real(dp) :: surf_reflect, b_surface, b_top, b_surface_SH4

    real(dp), dimension(nlay) :: beta, gam

    real(dp), dimension(nstr,nlay) :: eta, del
    real(dp), dimension(nlay) :: lam1, lam2, exptrm1, exptrm2
    real(dp), dimension(nlay) :: delta
    real(dp), dimension(nlay) ::  Q1, Q2, R1, R2, S1, S2
    real(dp), dimension(nlay) :: p1pl, p2pl, q1pl, q2pl, p1mn, p2mn, q1mn, q2mn
    real(dp), dimension(nlay) :: z1mn, z2mn, z1pl, z2pl
    real(dp), dimension(nlay) :: z1mn_up, z2mn_up, z1pl_up, z2pl_up, z1mn_down, z2mn_down, z1pl_down, z2pl_down
    real(dp), dimension(nlay) :: f00, f01, f02, f03, f10, f11, f12, f13, f20, f21, f22, f23, f30, f31, f32, f33
    real(dp), dimension(nlev) :: expon, tau_e

    real(dp), dimension(11,4*nlay) :: Mb
    real(dp), dimension(16,4*nlay) :: Mb_F
    real(dp), dimension(4*nlay) :: B

    real(dp), dimension(4*nlev, 4*nlay) :: F
    real(dp), dimension(4*nlev) :: G
    real(dp), dimension(4*nlay) :: F_bot
    real(dp) :: G_bot

    real(dp), dimension(4*nlev) :: flux_temp
    real(dp) :: flux_bot, f0
    real(dp), parameter :: eps_20 = 1.0e-20_dp


    real(dp), dimension(nlay) :: sigma_sq, pmom2, fc, c
    real(dp), dimension(nlay) :: w0, hg
    integer :: info
    integer, dimension(4*nlay) :: ipiv

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
    where (g_in(:) >= 1e-6_dp)
      fc(:) = g_in(:)**(nstr)
      pmom2(:) = g_in(:)**(nstr+1)
      sigma_sq(:) = real((nstr+1)**2 - nstr**2,dp) / &
      & ( log(fc(:)**2/pmom2(:)**2) )
      c(:) = exp(real(nstr**2,dp)/(2.0_dp*sigma_sq(:)))
      fc(:) = c(:)*fc(:)

      w0(:) = w_in(:)*((1.0_dp - fc(:))/(1.0_dp - fc(:)*w_in(:)))
      dtau(:) = (1.0_dp - w_in(:)*fc(:))*dtau(:)

    elsewhere
      w0(:) = w_in(:)
      fc(:) = 0.0_dp
    end where

    hg(:) = g_in(:)
    mu_z = mu_in(nlev) ! Can't do spherical geometry yet

    !! Reform edge optical depths
    tau_e(1) = 0.0_dp
    do k = 1, nlay
      tau_e(k+1) = tau_e(k) + dtau(k)
    end do

    f0 = 1.0_dp/mu_z

    Pu0(1) = 1.0_dp
    Pu0(2) = -mu_z
    Pu0(3) = (3.0_dp * mu_z**2 - 1.0_dp)/2.0_dp
    Pu0(4) = (5.0_dp * -mu_z**3 - 3.0_dp * -(mu_z))/2.0_dp

    w_multi(1,:) = 1.0_dp
    w_multi(2,:) = 3.0_dp * (hg(:) - fc(:)) / (1.0_dp - fc(:))
    w_multi(3,:) = 5.0_dp * (hg(:)**2 - fc(:)) / (1.0_dp - fc(:))
    w_multi(4,:) = 7.0_dp * (hg(:)**3 - fc(:)) / (1.0_dp - fc(:))

    do l = 1, nstr
      a(l,:) = real(2*(l-1) + 1,dp) -  w0(:) * w_multi(l,:) + eps_20
      bb(l,:) = ((w0(:) * w_multi(l,:)) * 1.0_dp * Pu0(l)) / (fourpi)
    end do

    surf_reflect = 0.0_dp

    b_surface = 0.0_dp + surf_reflect*mu_z*1.0_dp*exp(-tau_e(nlev)/mu_z)
    b_surface_SH4 = -(0.0_dp + surf_reflect*mu_z*1.0_dp*exp(-tau_e(nlev)/mu_z))/4.0_dp

    b_top = 0.0_dp

    !! Find beta and gamma
    beta(:) = a(1,:)*a(2,:) + (4.0_dp/9.0_dp)*a(1,:)*a(4,:) + (1.0_dp/9.0_dp)*a(3,:)*a(4,:)
    gam(:) = (a(1,:)*a(2,:)*a(3,:)*a(4,:))/9.0_dp

    ! Find k values - lambda in Rooney
    lam1(:) = sqrt((beta(:) + sqrt((beta(:)**2 - 4.0_dp*gam(:))))/2.0_dp)
    lam2(:) = sqrt((beta(:) - sqrt((beta(:)**2 - 4.0_dp*gam(:))))/2.0_dp)

    !! Find the delta values
    delta(:) = 9.0_dp*(f0**4 - beta(:)*f0**2 + gam(:))
    del(1,:) = (a(2,:)*bb(1,:) - bb(2,:)*f0)*(a(3,:)*a(4,:) - 9.0_dp*f0**2) &
     & + 2.0_dp*f0**2*(a(4,:)*bb(3,:) - 2.0_dp*a(4,:)*bb(1,:) - 3.0_dp*bb(4,:)*f0)
    del(2,:) = (a(1,:)*bb(2,:) - bb(1,:)*f0)*(a(3,:)*a(4,:) - 9.0_dp*f0**2) & 
      & - 2.0_dp*a(1,:)*f0*(a(4,:)*bb(3,:) - 3.0_dp*bb(4,:)*f0)
    del(3,:) = (a(4,:)*bb(3,:) - 3.0_dp*bb(4,:)*f0)*(a(1,:)*a(2,:) & 
      & - f0**2) - 2.0_dp*a(4,:)*f0*(a(1,:)*bb(2,:) - bb(1,:)*f0)
    del(4,:) = (a(3,:)*bb(4,:) - 3.0_dp*bb(3,:)*f0)*(a(1,:)*a(2,:) - f0**2) &
      &  + 2.0_dp*f0**2*(3.0_dp*a(1,:)*bb(2,:) - 2.0_dp*a(1,:)*bb(4,:) - 3.0_dp*bb(1,:)*f0)

    eta(1,:) = del(1,:)/delta(:)
    eta(2,:) = del(2,:)/delta(:)
    eta(3,:) = del(3,:)/delta(:)
    eta(4,:) = del(4,:)/delta(:)


    z1pl(:) = (eta(1,:)/2.0_dp + eta(2,:) + 5.0_dp*eta(3,:)/8.0_dp) * twopi
    z1mn(:) = (eta(1,:)/2.0_dp - eta(2,:) + 5.0_dp*eta(3,:)/8.0_dp) * twopi
    z2pl(:) = (-eta(1,:)/8.0_dp + 5.0_dp*eta(3,:)/8.0_dp + eta(4,:)) * twopi 
    z2mn(:) = (-eta(1,:)/8.0_dp + 5.0_dp*eta(3,:)/8.0_dp - eta(4,:)) * twopi

    !! Find e values
    exptrm1(:) = min(lam1(:)*dtau(:),35.0_dp)
    exptrm1(:) = exp(-exptrm1(:))
    exptrm2(:) = min(lam2(:)*dtau(:),35.0_dp)
    exptrm2(:) = exp(-exptrm2(:))  

    R1(:) = -a(1,:)/lam1(:); R2(:) = -a(1,:)/lam2(:)
    Q1(:) = 1.0_dp/2.0_dp * (a(1,:)*a(2,:)/(lam1(:)**2) - 1.0_dp)
    Q2(:) = 1.0_dp/2.0_dp * (a(1,:)*a(2,:)/(lam2(:)**2) - 1.0_dp)
    S1(:) = -3.0_dp/(2.0_dp*a(4,:)) * (a(1,:)*a(2,:)/lam1(:) - lam1(:));
    S2(:) = -3.0_dp/(2.0_dp*a(4,:)) * (a(1,:)*a(2,:)/lam2(:) - lam2(:))

    p1pl(:) = (1.0_dp/2.0_dp + R1(:) + 5.0_dp*Q1(:)/8.0_dp)  * twopi
    p2pl(:) = (1.0_dp/2.0_dp + R2(:) + 5.0_dp*Q2(:)/8.0_dp)  * twopi
    q1pl(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q1(:)/8.0_dp + S1(:)) * twopi
    q2pl(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q2(:)/8.0_dp + S2(:)) * twopi
    p1mn(:) = (1.0_dp/2.0_dp - R1(:) + 5.0_dp*Q1(:)/8.0_dp)  * twopi
    p2mn(:) = (1.0_dp/2.0_dp - R2(:) + 5.0_dp*Q2(:)/8.0_dp)  * twopi
    q1mn(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q1(:)/8.0_dp - S1(:)) * twopi
    q2mn(:) = (-1.0_dp/8.0_dp + 5.0_dp*Q2(:)/8.0_dp - S2(:)) * twopi

    f00(:)= p1mn(:)*exptrm1(:); f01(:) = p1pl(:)/exptrm1(:); f02(:) = p2mn(:)*exptrm2(:); f03(:) = p2pl(:)/exptrm2(:)
    f10(:) = q1mn(:)*exptrm1(:); f11(:) = q1pl(:)/exptrm1(:); f12(:) = q2mn(:)*exptrm2(:); f13(:) = q2pl(:)/exptrm2(:)
    f20(:) = p1pl(:)*exptrm1(:); f21(:) = p1mn(:)/exptrm1(:); f22(:) = p2pl(:)*exptrm2(:); f23(:) = p2mn(:)/exptrm2(:)
    f30(:) = q1pl(:)*exptrm1(:); f31(:) = q1mn(:)/exptrm1(:); f32(:) = q2pl(:)*exptrm2(:); f33(:) = q2mn(:)/exptrm2(:)

    expon(:) = exp(-tau_e(:)/mu_z)
    z1mn_up(:) = z1mn(:) * expon(2:nlev)
    z2mn_up(:) = z2mn(:) * expon(2:nlev)
    z1pl_up(:) = z1pl(:) * expon(2:nlev)
    z2pl_up(:) = z2pl(:) * expon(2:nlev)
    z1mn_down(:) = z1mn(:) * expon(1:nlay)
    z2mn_down(:) = z2mn(:) * expon(1:nlay)
    z1pl_down(:) = z1pl(:) * expon(1:nlay)
    z2pl_down(:) = z2pl(:) * expon(1:nlay)

    Mb(:,:) = 0.0_dp
    B(:) = 0.0_dp
    F_bot(:) = 0.0_dp
    G_bot = 0.0_dp
    F(:,:) = 0.0_dp
    G(:) = 0.0_dp

    !top boundary conditions
    Mb(6,1) = p1mn(1)
    Mb(6,2) = q1pl(1)
    Mb(5,2) = p1pl(1)
    Mb(5,3) = q2mn(1)
    Mb(4,3) = p2mn(1)
    Mb(4,4) = q2pl(1)
    Mb(3,4) = p2pl(1)
    Mb(7,1) = q1mn(1)

    B(1) = b_top - z1mn_down(1)
    B(2) = -b_top/4.0_dp - z2mn_down(1)

    ! bottom boundary conditions
    Mb(6,4*nlay-1) = f22(nlay) - surf_reflect*f02(nlay)
    Mb(6,4*nlay) = f33(nlay) - surf_reflect*f13(nlay)
    Mb(5,4*nlay) = f23(nlay) - surf_reflect*f03(nlay)
    Mb(7,4*nlay-2) = f21(nlay) - surf_reflect*f01(nlay)
    Mb(7,4*nlay-1) = f32(nlay) - surf_reflect*f12(nlay)
    Mb(8,4*nlay-3) = f20(nlay) - surf_reflect*f00(nlay)
    Mb(8,4*nlay-2) = f31(nlay) - surf_reflect*f11(nlay)
    Mb(9,4*nlay-3) = f30(nlay) - surf_reflect*f10(nlay)

    B(4*nlay-1) = b_surface - z1pl_up(nlay) + surf_reflect*z1mn_up(nlay)
    B(4*nlay) = b_surface_SH4 - z2pl_up(nlay) + surf_reflect*z2mn_up(nlay)

    !fill remaining rows of matrix
    do i = 1, nlay-1
      Mb(6,4*i - 1) = f02(i)
      Mb(6,4*i) = f13(i)
      Mb(6,4*i + 1) = -p1pl(i+1)
      Mb(6,4*i + 2) = -q1mn(i+1)
        
      Mb(5,4*i) = f03(i)
      Mb(5,4*i + 1) = -q1mn(i+1)
      Mb(5,4*i + 2) = -p1mn(i+1)
      Mb(5,4*i + 3) = -q2pl(i+1)
        
      Mb(4,4*i + 1) = -p1mn(i+1)
      Mb(4,4*i + 2) = -q1pl(i+1)
      Mb(4,4*i + 3) = -p2pl(i+1)
      Mb(4,4*i + 4) = -q2mn(i+1)
        
      Mb(3,4*i + 2) = -p1pl(i+1)
      Mb(3,4*i + 3) = -q2mn(i+1)
      Mb(3,4*i + 4) = -p2mn(i+1)
        
      Mb(2,4*i + 3) = -p2mn(i+1)
      Mb(2,4*i + 4) = -q2pl(i+1)
        
      Mb(1,4*i + 4) = -p2pl(i+1)
        
      Mb(7,4*i - 2) = f01(i)
      Mb(7,4*i -1) = f12(i)
      Mb(7,4*i) = f23(i)
      Mb(7,4*i + 1) = -q1pl(i+1)
        
      Mb(8,4*i - 3) = f00(i)
      Mb(8,4*i - 2) = f11(i)
      Mb(8,4*i - 1) = f22(i)
      Mb(8,4*i) = f33(i)
        
      Mb(9,4*i - 3) = f10(i)
      Mb(9,4*i - 2) = f21(i)
      Mb(9,4*i - 1) = f32(i)
        
      Mb(10,4*i - 3) = f20(i)
      Mb(10,4*i - 2) = f31(i)
        
      Mb(11,4*i - 3) = f30(i)
    end do 

    do i = 1, nlay-1
      B(4*i - 1) = z1mn_down(i+1) - z1mn_up(i)
      B(4*i) = z2mn_down(i+1) - z2mn_up(i)
      B(4*i + 1) = z1pl_down(i+1)- z1pl_up(i)
      B(4*i + 2) = z2pl_down(i+1) - z2pl_up(i)
    end do

    Mb_F(:,:) = 0.0_dp
    Mb_F(6,:) = Mb(1,:)
    Mb_F(7,:) = Mb(2,:)
    Mb_F(8,:) = Mb(3,:)  
    Mb_F(9,:) = Mb(4,:)
    Mb_F(10,:) = Mb(5,:)  
    Mb_F(11,:) = Mb(6,:)
    Mb_F(12,:) = Mb(7,:)
    Mb_F(13,:) = Mb(8,:)
    Mb_F(14,:) = Mb(9,:)   
    Mb_F(15,:) = Mb(10,:)    
    Mb_F(16,:) = Mb(11,:)

    call dgbsv(4*nlay, 5, 5, 1, Mb_F, 16, ipiv, B, 4*nlay, info)

    if (info == 0) then
        ! Success, B now contains the solution
    else
        print *, "An error occurred: ", info
    endif

   ! flux at bottom of atmosphere
    F_bot(nlay-3) = f20(nlay)
    F_bot(nlay-2) = f21(nlay)
    F_bot(nlay-1) = f22(nlay)
    F_bot(nlay) = f23(nlay)
    G_bot = z1pl_up(nlay)

    F(1,1) = p1mn(1)
    F(1,2) = p1pl(1)
    F(1,3) = p2mn(1)
    F(1,4) = p2pl(1)
    F(2,1) = q1mn(1)
    F(2,2) = q1pl(1)
    F(2,3) = q2mn(1)
    F(2,4) = q2pl(1)
    F(3,1) = p1pl(1)
    F(3,2) = p1mn(1)
    F(3,3) = p2pl(1)
    F(3,4) = p2mn(1)
    F(4,1) = q1pl(1)
    F(4,2) = q1mn(1)
    F(4,2) = q2pl(1)
    F(4,4) = q2mn(1)

    k = 0
    do i = 1, 4*nlay, 4
      F(i+4,i) = f00(k+1)
      F(i+4,i+1) = f01(k+1)
      F(i+4,i+2) = f02(k+1)
      F(i+4,i+3) = f03(k+1)
      F(i+5,i) = f10(k+1)
      F(i+5,i+1) = f11(k+1)
      F(i+5,i+2) = f12(k+1)
      F(i+5,i+3) = f13(k+1)
      F(i+6,i) = f20(k+1)
      F(i+6,i+1) = f21(k+1)
      F(i+6,i+2) = f22(k+1)
      F(i+6,i+3) = f23(k+1)
      F(i+7,i) = f30(k+1)
      F(i+7,i+1) = f31(k+1)
      F(i+7,i+2) = f32(k+1)
      F(i+7,i+3) = f33(k+1)
      k = k + 1
    end do
                             
    G(1) = z1mn_down(1)
    G(2) = z2mn_down(1)
    G(3) = z1pl_down(1)
    G(4) = z2pl_down(1)
    do i = 1, nlay
      G(4*i + 1) = z1mn_up(i)
      G(4*i + 2) = z2mn_up(i)
      G(4*i + 3) = z1pl_up(i)
      G(4*i + 4) = z2pl_up(i)
    end do

    flux_temp(:) = matmul(F(:,:),B(:)) + G(:)
    flux_bot = sum(F_bot(:)*B(:)) + G_bot

    do i = 1, nlev
      flx_down(i) = max(flux_temp(i*4 - 3),0.0_dp)*F0_in + mu_z * F0_in * expon(i)
      flx_up(i) = max(flux_temp(i*4 - 1), 0.0_dp)*F0_in
    end do

  end subroutine sw_SH_four_stream

end module sw_SH4_mod
