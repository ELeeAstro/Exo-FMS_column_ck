!!!
! Elspeth KH Lee - Jun 2024 : Overhauled
!
! sw: Adding method for shortwave radiation, used in various Earth GCMs (e.g. Mendonça et al. 2015)
!     Pros: Very fast method with approximate scattering, decently accurate
!     Cons: Not multiple scattering
!!!


module sw_adding_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  private :: sw_approx_adding
  public :: sw_adding

contains

  subroutine sw_adding(nlay, nlev, nb, ng, gw, tau_e, mu_z, Finc, ssa, gg, a_surf, sw_up, sw_down, sw_net, asr)
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
          call sw_approx_adding(nlay, nlev, Finc(b), mu_z(:), tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), a_surf(b), & 
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

  end subroutine sw_adding


  subroutine sw_approx_adding(nlay, nlev, F0_in, mu_in, tau_in, w_in, g_in, w_surf_in, flx_down, flx_up)
    implicit none

    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0_in, w_surf_in
    real(dp), dimension(nlev), intent(in) :: tau_in, mu_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    integer :: k
    real(dp) :: lamtau, e_lamtau, arg, apg, amg
    real(dp), dimension(nlev) ::  w, g, f
    real(dp), dimension(nlev) :: tau_e
    real(dp), dimension(nlay) :: tau
    real(dp), dimension(nlev) :: tau_s, w_s, g_s
    real(dp), dimension(nlev) :: lam, u, N, gam, alp
    real(dp), dimension(nlev) :: R_b, T_b, R, T
    real(dp), dimension(nlev) :: Tf
    real(dp), dimension(nlev) :: cum_trans

    !! Design w and g to include surface property level
    w(1:nlay) = w_in(:)
    g(1:nlay) = g_in(:)

    w(nlev) = w_surf_in
    g(nlev) = 0.0_dp

    ! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w(:) <= 1.0e-12_dp)) then

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

    !! Backscattering approximation
    f(:) = g(:)**2

    !! Do optical depth rescaling
    tau_e(1) = tau_in(1)
    do k = 1, nlay
      tau(k) = tau_in(k+1) - tau_in(k)
      tau_s(k) = tau(k) * (1.0_dp - w(k)*f(k))
      tau_e(k+1) = tau_e(k) + tau_s(k)
    end do

    !! Perform approximate adding method
    do k = 1, nlev

      w_s(k) = w(k) * ((1.0_dp - f(k))/(1.0_dp - w(k)*f(k)))
      g_s(k) = (g(k) - f(k))/(1.0_dp - f(k))
      lam(k) = sqrt(3.0_dp*(1.0_dp - w_s(k))*(1.0_dp - w_s(k)*g_s(k)))
      gam(k) = 0.5_dp * w_s(k) * (1.0_dp + 3.0_dp*g_s(k)*(1.0_dp - w_s(k))*mu_in(k)**2)/(1.0_dp - lam(k)**2*mu_in(k)**2)
      alp(k) = 0.75_dp * w_s(k) * mu_in(k) * (1.0_dp + g_s(k)*(1.0_dp - w_s(k)))/(1.0_dp - lam(k)**2*mu_in(k)**2)
      u(k) = (3.0_dp/2.0_dp) * ((1.0_dp - w_s(k)*g_s(k))/lam(k))

      lamtau = min(lam(k)*tau_e(k),99.0_dp)
      e_lamtau = exp(-lamtau)

      N(k) = (u(k) + 1.0_dp)**2 * 1.0_dp/e_lamtau - (u(k) - 1.0_dp)**2  * e_lamtau

      R_b(k) = (u(k) + 1.0_dp)*(u(k) - 1.0_dp)*(1.0_dp/e_lamtau - e_lamtau)/N(k)
      T_b(k) = 4.0_dp * u(k)/N(k)

      arg = min(tau_e(k)/mu_in(k),99.0_dp)
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
      flx_down(k) = Tf(k) + ((T(k) - Tf(k)) +  &
      & Tf(k)*R(k+1)*R_b(k))/(1.0_dp - R_b(k)*R_b(k+1))
    end do
    flx_down(nlev) = Tf(nlev)

    !! Calculate upward flux
    do k = 1, nlay
      flx_up(k) = (Tf(k)*R(k+1) + (T(k) - Tf(k))*R_b(k+1))/(1.0_dp - R_b(k)*R_b(k+1))
    end do
    flx_up(nlev) = flx_down(nlev) * w(nlev)

    !! Scale with the incident flux
    flx_down(:) = flx_down(:) * mu_in(nlev) * F0_in
    flx_up(:) = flx_up(:) * mu_in(nlev) * F0_in

  end subroutine sw_approx_adding

end module sw_adding_mod
