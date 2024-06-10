!!!
! Elspeth KH Lee - May 2021 : Initial version
!                - Jan 2022 : Working version
!                - Jun 2024 : Optimisations
!
! lw: Two-stream method following the "Toon89" method (Toon et al. 1989)
!     Based on the CHIMERA code by Mike Line, but cleaned up slightly
!     Pros: Fast, accurate at high optical depths, well used and familiar method
!     Cons: For longwave combined high ssa and g (0.9+,0.9+) can be unstable
!!!

module lw_Toon_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4  
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp
  real(dp), parameter :: hp = 6.62607015e-34_dp
  real(dp), parameter :: kb = 1.380649e-23_dp
  real(dp), parameter :: c_s = 2.99792458e8_dp
  real(dp), parameter :: c1 = (hp * c_s) / kb
  real(dp), parameter :: c2 = c_s**2
  real(dp), parameter :: n2 = 2.0_dp * hp * c2

  real(dp), parameter :: ubari = 0.5_dp

  !! Gauss quadrature variables, cosine angle values (uarr) and weights (w)
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these
  
  !! Optimised quadrature for 1 node (Hogan 2024)
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.6096748751_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! Legendre quadrature for 2 nodes
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.21132487_dp, 0.78867513_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  private :: lw_Toon89, dtridgl, BB_integrate
  public :: lw_Toon

contains


  subroutine lw_Toon(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, a_surf, Tint, lw_up, lw_down, lw_net, olr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev, nb, ng
    real(dp), dimension(ng), intent(in) :: gw
    real(dp), dimension(nb+1), intent(in) :: wn_e
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(ng,nb,nlev), intent(in) :: tau_e
    real(dp), dimension(ng,nb,nlay), intent(in) :: ssa, gg
    real(dp), dimension(nb), intent(in) :: a_surf
    real(dp), intent(in) :: Tint

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down, lw_net
    real(dp), intent(out) :: olr

    !! Work variables
    integer :: k, b, g
    real(dp), dimension(nlev) :: Te
    real(dp), dimension(nb,nlev) :: be
    real(dp), dimension(nb) :: be_int
    real(dp), dimension(nb,nlev) :: lw_down_b, lw_up_b
    real(dp), dimension(ng,nlev) :: lw_down_g, lw_up_g

    !! Use WENO4 method to (smoothly) interpolate layers to levels
    Te(:) = interpolate_weno4(pe, pl, Tl, .False.)

    !! Edges are linearly interpolated to avoid overshoot
    Te(1) = 10.0_dp**(log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Tl(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    !! Find integrated planck function for each band for each level
    do k = 1, nlev
      call BB_integrate(nb, Te(k), wn_e, be(:,k))   ! Integrated planck function intensity at levels for each band
    end do
    call BB_integrate(nb, Tint, wn_e, be_int) ! Integrated planck function intensity for internal temperature

    !! Longwave flux calculation
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp
    do b = 1, nb
      lw_up_b(b,:) = 0.0_dp
      lw_down_b(b,:) = 0.0_dp
      do g = 1, ng
        call lw_Toon89(nlay, nlev, be(b,:), tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), a_surf(b), &
          & lw_up_g(g,:), lw_down_g(g,:))
        lw_up_b(b,:) = lw_up_b(b,:) + lw_up_g(g,:) * gw(g)
        lw_down_b(b,:) = lw_down_b(b,:) + lw_down_g(g,:) * gw(g)
      end do
      lw_up(:) = lw_up(:) + lw_up_b(b,:)
      lw_down(:) = lw_down(:) + lw_down_b(b,:)
    end do

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Apply net internal flux at lower boundary
    lw_net(nlev) = sb * Tint**4

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_Toon

  subroutine lw_Toon89(nlay, nlev, be, tau_in, w_in, g_in, a_surf_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: a_surf_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables
    integer :: k, i, n, m
    integer :: l, lm2, lm1
    real(dp) :: Bsurf, Btop, bottom, tautop
    real(dp), dimension(nlev) :: tau
    real(dp), dimension(nlay) :: dtau_in, dtau
    real(dp), dimension(nlay) :: w0, hg
    real(dp), dimension(nlay) :: B0, B1
    real(dp), dimension(nlay) :: lam, gam, alp, term
    real(dp), dimension(nlay) :: Cpm1, Cmm1, Cp, Cm
    real(dp), dimension(nlay) :: exptrm, Ep, Em, E1, E2, E3, E4
    real(dp), dimension(nlay+nlay) :: Af, Bf, Cf, Df, xkk
    real(dp), dimension(nlay) :: xk1, xk2

    real(dp), dimension(nlay) :: g, h, xj, xk
    real(dp), dimension(nlay) :: alpha1, alpha2, sigma1, sigma2
    real(dp), dimension(nlay) :: em1, obj, epp, obj2, epp2, em2, em3

    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    l = nlay + nlay
    lm2 = l - 2
    lm1 = l - 1

    dtau_in(:) = tau_in(2:nlev) - tau_in(1:nlay)

    ! Delta eddington scaling
    w0(:) = (1.0_dp - g_in(:)**2)*w_in(:)/(1.0_dp - w_in(:)*g_in(:)**2)
    dtau(:) = (1.0_dp - w_in(:)*g_in(:)**2)*dtau_in(:)
    hg(:) = g_in(:)/(1.0_dp + g_in(:))

    tau(1) = 0.0_dp
    do k = 1, nlay
      tau(k+1) = tau(k) + dtau(k)
    end do

    alp(:) = sqrt((1.0_dp - w0(:))/(1.0_dp - w0(:)*hg(:)))
    lam(:) = alp(:)*(1.0_dp - w0(:)*hg(:))/ubari
    gam(:) = (1.0_dp - alp(:))/(1.0_dp + alp(:))
    term(:) = ubari/(1.0_dp - w0(:)*hg(:))

    do k = 1, nlay
      if (dtau(k) <= 1.0e-6_dp) then
        ! For low optical depths use the isothermal approimation
        B1(k) = 0.0_dp
        B0(k) = 0.5_dp*(be(k+1) + be(k))
      else
        B1(k) = (be(k+1) - be(k))/dtau(k) ! Linear in tau term
        B0(k) = be(k)
      endif
    end do

    !Cpm1 and Cmm1 are the C+ and C- terms evaluated at the top of the layer.
    Cpm1(:) = B0(:) + B1(:)*term(:)
    Cmm1(:) = B0(:) - B1(:)*term(:)
    !Cp and Cm are the C+ and C- terms evaluated at the bottom of the layer.
    Cp(:) = B0(:) + B1(:)*dtau(:) + B1(:)*term(:)
    Cm(:) = B0(:) + B1(:)*dtau(:) - B1(:)*term(:)

    tautop = dtau(1)*exp(-1.0_dp)
    Btop = (1.0_dp - exp(-tautop/ubari))*be(1)
    Bsurf = be(nlev)
    bottom = Bsurf + B1(nlay)*ubari

    !Solve for the coefficients of system of equations using boundary conditions
    !Exponential terms:
    exptrm(:) = min(lam(:)*dtau(:),35.0_dp)
    Ep(:) = exp(exptrm(:))
    Em(:) = 1.0_dp/Ep(:)

    E1(:) = Ep(:) + gam(:)*Em(:)
    E2(:) = Ep(:) - gam(:)*Em(:)
    E3(:) = gam(:)*Ep(:) + Em(:)
    E4(:) = gam(:)*Ep(:) - Em(:)

    Af(1) = 0.0_dp
    Bf(1) = gam(1) + 1.0_dp
    Cf(1) = gam(1) - 1.0_dp
    Df(1) = btop - Cmm1(1)

    n = 0
    do i = 2, lm2, 2
      n = n + 1
      Af(i) = (E1(n)+E3(n))*(gam(n+1) - 1.0_dp)
      Bf(i) = (E2(n)+E4(n))*(gam(n+1) - 1.0_dp)
      Cf(i) = 2.0_dp*(1.0_dp - gam(n+1)**2)
      Df(i) = (gam(n+1) - 1.0_dp)*(Cpm1(n+1) - Cp(n)) + (1.0_dp - gam(n+1))*(Cm(n) - Cmm1(n+1))
    end do

    n = 0
    do i = 3, lm1, 2
      n = n + 1
      Af(i) = 2.0_dp*(1.0_dp - gam(n)**2)
      Bf(i) = (E1(n)-E3(n))*(1.0_dp + gam(n+1))
      Cf(i) = (E1(n)+E3(n))*(gam(n+1)-1.0_dp)
      Df(i) = E3(n)*(Cpm1(n+1) - Cp(n)) + E1(n)*(Cm(n) - Cmm1(n+1))
    end do

    Af(l) = E1(nlay) - a_surf_in*E3(nlay)
    Bf(l) = E2(nlay) - a_surf_in*E4(nlay)
    Cf(l) = 0.0_dp
    Df(l) = bsurf - Cp(nlay) + a_surf_in*Cm(nlay)

    call dtridgl(l, Af, Bf, Cf, Df, xkk)

    do n = 1, nlay
      xk1(n) = xkk(2*n-1) + xkk(2*n)
      xk2(n) = xkk(2*n-1) - xkk(2*n)
      if (abs(xk2(n)/xkk(2*n-1)) < 1e-30_dp) then
        xk2(n) = 0.0_dp
      end if
      !print*, xk1(n), xk2(n)
    end do

    where (w0(:) <= 1e-4_dp)
      g(:) = 0.0_dp
      h(:) = 0.0_dp
      xj(:) = 0.0_dp
      xk(:) = 0.0_dp
      alpha1(:)=twopi*B0(:)
      alpha2(:)=twopi*B1(:)
      sigma1(:)=alpha1(:)
      sigma2(:)=alpha2(:)
    elsewhere 
      g(:)=twopi*w0(:)*xk1(:)*(1.0_dp+hg(:)*alp(:))/(1.0_dp+alp(:))
      h(:)=twopi*w0(:)*xk2(:)*(1.0_dp-hg(:)*alp(:))/(1.0_dp+alp(:))
      xj(:)=twopi*w0(:)*xk1(:)*(1.0_dp-hg(:)*alp(:))/(1.0_dp+alp(:))
      xk(:)=twopi*w0(:)*xk2(:)*(1.0_dp+hg(:)*alp(:))/(1.0_dp+alp(:))
      alpha1(:)=twopi*(B0(:)+B1(:)*(ubari*w0(:)*hg(:)/(1.0_dp-w0(:)*hg(:))))
      alpha2(:)=twopi*B1(:)
      sigma1(:)=twopi*(B0(:)-B1(:)*(ubari*w0(:)*hg(:)/(1.0_dp-w0(:)*hg(:))))
      sigma2(:)=alpha2(:)
    end where

    obj(:) = min(lam(:)*dtau(:),35.0_dp)
    epp(:) = exp(obj(:))
    em1(:) = 1.0_dp/epp(:)
    obj2(:) = min(0.5_dp*lam(:)*dtau(:),35.0_dp)
    epp2(:) = exp(obj2(:))

    ! Zero the total flux arrays
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp

    do m = 1, nmu

      em2(:) = exp(-dtau(:)/uarr(m))
      em3(:) = em1(:)*em2(:)

      !! Begin two-stream loops
      !! Perform downward loop first
      ! Top boundary condition - intensity downward from top boundary (tautop, assumed isothermal)
      lw_down_g(1) = twopi*(1.0_dp - exp(-tautop/uarr(m)))*be(1)
      do k = 1, nlay
        lw_down_g(k+1) = lw_down_g(k)*em2(k) + &
        & (xj(k)/(lam(k)*uarr(m)+1.0_dp))*(epp(k)-em2(k)) + &
        & (xk(k)/(lam(k)*uarr(m)-1.0_dp))*(em2(k)-em(k))+sigma1(k)*(1.0_dp-em2(k)) + &
        & sigma2(k)*(uarr(m)*em2(k)+dtau(k)-uarr(m))
      end do

      lw_up_g(nlev) = twopi*(Bsurf+B1(nlay)*uarr(m))
      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*em2(k) + &
        & (g(k)/(lam(k)*uarr(m)-1.0_dp))*(epp(k)*em2(k)-1.0_dp) + &
        & (h(k)/(lam(k)*uarr(m)+1.0_dp))*(1.0_dp-em3(k))+alpha1(k)*(1.0_dp-em2(k)) + &
        & alpha2(k)*(uarr(m)-(dtau(k)+uarr(m))*em2(k))
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      flx_down(:) = flx_down(:) + lw_down_g(:) * wuarr(m)
      flx_up(:) = flx_up(:) + lw_up_g(:) * wuarr(m)

    end do

  end subroutine lw_Toon89

  subroutine dtridgl(l, af, bf, cf, df, xk)
    implicit none

    integer, intent(in) :: l
    real(dp), dimension(l), intent(in) :: af, bf, cf, df
    real(dp), dimension(l), intent(out) :: xk

    integer :: i
    integer, parameter :: nmax = 301
    real(dp) :: x
    real(dp), dimension(nmax) :: as, ds

    as(l) = af(l)/bf(l)
    ds(l) = df(l)/bf(l)

    do i = l-1, 1, -1
      x = 1.0_dp/(bf(i) - cf(i)*as(i+1))
      as(i) = af(i)*x
      ds(i) = (df(i) - cf(i)*ds(i+1))*x
    end do

    xk(1) = ds(1)

    do i = 2, l
      xk(i) = ds(i)-as(i)*xk(i-1)
    end do

  end subroutine dtridgl

  subroutine BB_integrate(n_b, Te, wn_e, be)
    implicit none

    integer, intent(in) :: n_b
    real(dp), intent(in) :: Te
    real(dp), dimension(n_b+1), intent(in) :: wn_e

    real(dp), dimension(n_b), intent(out) :: be

    integer :: ww, j, intitera
    real(dp), dimension(n_b+1) :: iB
    real(dp) :: x, x2, x3, itera, summ, dn

    !! Code for integrating the blckbody function between two wavenumbers
    !! This is a method that uses a sum convergence
    !! Taken from: spectralcalc.com/blackbody/inband_radiance.html

      if (Te < 1e-6_dp) then
        be(:) = 0.0_dp
        return
      end if

      do ww = 1, n_b+1

        x = c1 * 100.0_dp * wn_e(ww)/ Te
        x2 = x**2
        x3 = x**3

        itera = 2.0_dp + 20.0_dp/x
        if (itera > 150) then
          itera = 150
        end if
        intitera = int(itera)

        summ = 0.0_dp
        do j = 1, intitera + 1
          dn = 1.0_dp/real(j,kind=dp)
          summ = summ +  exp(-min(real(j,kind=dp)*x,300.0_dp))* &
          & (x3 + (3.0_dp * x2 + 6.0_dp*(x+dn)*dn)*dn)*dn
        end do

        iB(ww) = n2 * (Te/c1)**(4) * summ
      end do

      do ww = 1, n_b
        be(ww) = max(iB(ww+1) - iB(ww),0.0_dp)
      end do

  end subroutine BB_integrate

end module lw_Toon_mod