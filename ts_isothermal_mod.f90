!!!
! Elspeth KH Lee - May 2021
! Two-stream method following the isothermal layer approximation
! Pros: Very fast
! Cons: Inaccurate at high optical depths (has no linear in tau term)
!!!

module ts_isothermal_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Common constants
  real(kind=dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
  real(kind=dp), parameter :: twopi = 2.0_dp*pi, fourpi = 4.0_dp*pi
  real(kind=dp), parameter :: sb = 5.670374419e-8_dp
  real(kind=dp), parameter :: h = 6.62607015e-34_dp
  real(kind=dp), parameter :: kb = 1.380649e-23_dp
  real(kind=dp), parameter :: c_s = 2.99792458e8_dp
  real(kind=dp), parameter :: c1 = (h * c_s) / kb
  real(kind=dp), parameter :: c2 = c_s**2
  real(kind=dp), parameter :: n2 = 2.0_dp * h * c2

  real(dp), parameter :: D = 1.66_dp  ! Diffusivity factor

  public :: ts_isothermal
  private :: lw_grey_updown, sw_grey_updown_adding

contains

  subroutine ts_isothermal(nlay, nlev, nb, ng, gw, wn_e, Tl, tau_e, ssa, gg, mu_z, Finc, Tint, olr, net_F)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev, nb, ng
    real(dp), dimension(ng), intent(in) :: gw
    real(dp), dimension(nlay), intent(in) :: Tl
    real(dp), dimension(ng,nb,nlev), intent(in) :: tau_e
    real(dp), dimension(ng,nb,nlay), intent(in) :: ssa, gg
    real(dp), dimension(nb), intent(in) :: Finc
    real(dp), dimension(nb+1), intent(in) :: wn_e
    real(dp), intent(in) :: mu_z, Tint

    !! Output variables
    real(dp), intent(out) :: olr
    real(dp), dimension(nlev), intent(out) :: net_F

    !! Work variables
    integer :: i, b, g
    real(dp), dimension(nb) :: be_int
    real(dp), dimension(nb,nlay) :: bl
    real(dp), dimension(nlev) :: sw_down, sw_up, lw_down, lw_up
    real(dp), dimension(nb,nlev) :: sw_up_b, sw_down_b
    real(dp), dimension(ng,nlev) :: sw_up_g, sw_down_g
    real(dp), dimension(nb,nlev) :: lw_down_b, lw_up_b
    real(dp), dimension(nlev) :: lw_net, sw_net

    real(dp) :: a_surf = 0.0_dp

    !! Shortwave flux calculation
    if (mu_z > 0.0_dp) then
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
      do b = 1, nb
        sw_down_b(b,:) = 0.0_dp
        sw_up_b(b,:) = 0.0_dp
        do g = 1, ng
          call sw_grey_updown_adding(nlay, nlev, Finc(b), tau_e(g,b,:), mu_z, ssa(g,b,:), gg(g,b,:), a_surf, &
          & sw_down_g(g,:), sw_up_g(g,:))
          sw_down_b(b,:) = sw_down_b(b,:) + sw_down_g(g,:) * gw(g)
          sw_up_b(b,:) = sw_up_b(b,:) + sw_up_g(g,:) * gw(g)
        end do
        !print*, b, sw_down_b(b,:)
        sw_down(:) = sw_down(:) + sw_down_b(b,:)
        sw_up(:) = sw_up(:) + sw_up_b(b,:)
      end do
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Longwave two-stream flux calculation
    do i = 1, nlay
      call BB_integrate(nlay, nb, Tl(i), wn_e(:), bl(:,i))
    end do
    call BB_integrate(nlev, nb, Tint, wn_e(:), be_int(:))

    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp
    do b = 1, nb
      call lw_grey_updown(nlay, nlev, ng, gw(:), bl(b,:), be_int(b), tau_e(:,b,:), lw_up_b(b,:), lw_down_b(b,:))
      lw_up(:) = lw_up(:) + lw_up_b(b,:)
      lw_down(:) = lw_down(:) + lw_down_b(b,:)
    end do

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

    !! Output olr
    olr = lw_up(1)

  end subroutine ts_isothermal

  subroutine lw_grey_updown(nlay, nlev, ng, gw, bl, be_int, tau_IRe, lw_up, lw_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev, ng
    real(dp), dimension(ng), intent(in) :: gw
    real(dp), dimension(nlay), intent(in) :: bl
    real(dp), dimension(ng,nlev), intent(in) :: tau_IRe
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k, g
    real(dp), dimension(nlay) :: dtau, Tp, B0
    real(dp), dimension(ng,nlay) :: lw_up_g, lw_down_g


    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp
    do g = 1, ng

      !! Prepare loop
      do k = 1, nlay
        dtau(k) = tau_IRe(g,k+1) - tau_IRe(g,k)
        Tp(k) = exp(-D*dtau(k))
      end do

      !! First do the downward loop
      lw_down_g(g,1) = 0.0_dp
      do k = 1, nlay
         lw_down_g(g,k+1) = lw_down_g(g,k)*Tp(k) + bl(k)*(1.0_dp - Tp(k))
      end do

      !! Perform upward loop
      ! Lower boundary condition - internal heat definition Fint = F_up - F_down
      lw_up_g(g,nlev) = lw_down_g(g,nlev) + be_int
      do k = nlay, 1, -1
        lw_up_g(g,k) = lw_up_g(g,k+1)*Tp(k) + bl(k)*(1.0_dp - Tp(k))
      end do

      lw_up(:) = lw_up(:) + lw_up_g(g,:) * gw(g)
      lw_down(:) = lw_down(:) + lw_down_g(g,:) * gw(g)

    end do


  end subroutine lw_grey_updown
  
  subroutine sw_grey_updown_adding(nlay, nlev, Finc, tau_Ve, mu_z, w_in, g_in, w_surf, sw_down, sw_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: Finc, mu_z, w_surf
    real(dp), dimension(nlev), intent(in) :: tau_Ve
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_down, sw_up

    !! Work variables
    integer :: k
    real(dp) :: lamtau, e_lamtau, lim, arg, apg, amg
    real(dp), dimension(nlev) ::  w, g, f
    real(dp), dimension(nlev) :: tau_Ve_s
    real(dp), dimension(nlay) :: tau
    real(dp), dimension(nlev) :: tau_s, w_s, f_s, g_s
    real(dp), dimension(nlev) :: lam, u, N, gam, alp
    real(dp), dimension(nlev) :: R_b, T_b, R, T
    real(dp), dimension(nlev) :: Tf

    ! Design w and g to include surface property level
    w(1:nlay) = w_in(:)
    g(1:nlay) = g_in(:)

    w(nlev) = 0.0_dp
    g(nlev) = 0.0_dp

    ! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w(:) <= 1.0e-12_dp)) then
      sw_down(:) = Finc * mu_z * exp(-tau_Ve(:)/mu_z)
      sw_down(nlev) = sw_down(nlev) * (1.0_dp - w_surf) ! The surface flux for surface heating is the amount of flux absorbed by surface
      sw_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo
      return
    end if

    w(nlev) = w_surf
    g(nlev) = 0.0_dp

    ! Backscattering approximation
    f(:) = g(:)**2

       !! Do optical depth rescaling
    tau_Ve_s(1) = tau_Ve(1)
    do k = 1, nlay
      tau(k) = tau_Ve(k+1) - tau_Ve(k)
      tau_s(k) = tau(k) * (1.0_dp - w(k)*f(k))
      tau_Ve_s(k+1) = tau_Ve_s(k) + tau_s(k)
    end do

    do k = 1, nlev

      w_s(k) = w(k) * ((1.0_dp - f(k))/(1.0_dp - w(k)*f(k)))
      g_s(k) = (g(k) - f(k))/(1.0_dp - f(k))
      lam(k) = sqrt(3.0_dp*(1.0_dp - w_s(k))*(1.0_dp - w_s(k)*g_s(k)))
      gam(k) = 0.5_dp * w_s(k) * (1.0_dp + 3.0_dp*g_s(k)*(1.0_dp - w_s(k))*mu_z**2)/(1.0_dp - lam(k)**2*mu_z**2)
      alp(k) = 0.75_dp * w_s(k) * mu_z * (1.0_dp + g_s(k)*(1.0_dp - w_s(k)))/(1.0_dp - lam(k)**2*mu_z**2)
      u(k) = (3.0_dp/2.0_dp) * ((1.0_dp - w_s(k)*g_s(k))/lam(k))

      lamtau = min(lam(k)*tau_Ve_s(k),99.0_dp)
      e_lamtau = exp(-lamtau)

      N(k) = (u(k) + 1.0_dp)**2 * 1.0_dp/e_lamtau - (u(k) - 1.0_dp)**2  * e_lamtau

      R_b(k) = (u(k) + 1.0_dp)*(u(k) - 1.0_dp)*(1.0_dp/e_lamtau - e_lamtau)/N(k)
      T_b(k) = 4.0_dp * u(k)/N(k)

      arg = min(tau_Ve_s(k)/mu_z,99.0_dp)
      Tf(k) = exp(-arg)

      apg = alp(k) + gam(k)
      amg = alp(k) - gam(k)

      R(k) = amg*(T_b(k)*Tf(k) - 1.0_dp) + apg*R_b(k)

      T(k) = apg*T_b(k) + (amg*R_b(k) - (apg - 1.0_dp))*Tf(k)

      R(k) = max(R(k), 0.0_dp)
      T(k) = max(T(k), 0.0_dp)
      R_b(k) = max(R_b(k), 0.0_dp)
      T_b(k) = max(T_b(k), 0.0_dp)

    end do

    !! Calculate downward flux
    do k = 1, nlay
      sw_down(k) = Tf(k) + ((T(k) - Tf(k)) +  &
      & Tf(k)*R(k+1)*R_b(k))/(1.0_dp - R_b(k)*R_b(k+1))
    end do
    sw_down(nlev) = Tf(nlev)

    !! Calculate upward flux
    do k = 1, nlay
      sw_up(k) = (Tf(k)*R(k+1) + (T(k) - Tf(k))*R_b(k+1))/(1.0_dp - R_b(k)*R_b(k+1))
    end do
    sw_up(nlev) = sw_down(nlev) * w_surf

    !! Scale with the incident flux
    sw_down(:) = sw_down(:) * mu_z * Finc
    sw_up(:) = sw_up(:) * mu_z * Finc

  end subroutine sw_grey_updown_adding
   
  subroutine BB_integrate(nlev, nb, Te, wn_e, be)
    implicit none

    integer, intent(in) :: nlev, nb
    real(dp), intent(in) :: Te
    real(dp), dimension(nb+1), intent(in) :: wn_e
    real(dp), dimension(nb), intent(out) :: be

    integer :: i, b, j, intitera
    real(dp), dimension(nb+1) :: iB
    real(dp) :: x, x2, x3, itera, summ, dn

    !! Code for integrating the blckbody function between two wavenumbers
    !! This is a method that uses a sum convergence
    !! Taken from: spectralcalc.com/blackbody/inband_radiance.html

    if (Te < 1.0e-6_dp) then

      be(:) = 1.0e-6_dp

    else

      do b = nb+1, 1, -1

        x = c1 * 100.0_dp * wn_e(b)/ Te
        x2 = x**2
        x3 = x**3

        itera = 2.0_dp + 20.0_dp/x
        if (itera > 512) then
          itera = 512
        end if
        intitera = int(itera)

        summ = 0.0_dp
        do j = 1, intitera + 1
          dn = 1.0_dp/real(j,dp)
          summ = summ + exp(-min(real(j,dp)*x,300.0_dp)) * &
          & (x3 + (3.0_dp * x2 + 6.0_dp*(x+dn)*dn)*dn)*dn
        end do

        iB(b) = n2 * (Te/c1)**4 * summ
        
      end do

      do b = nb, 1, -1
        be(b) = max(iB(b+1) - iB(b), 0.0_dp)
        !print*, b, be(b), itera, Te, wn_e(b)
      end do

    end if

  end subroutine BB_integrate

end module ts_isothermal_mod
