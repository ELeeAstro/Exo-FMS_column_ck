!!!
! Elspeth KH Lee - May 2021
! Two-stream method following the Toon method (Toon et al. 1989)
!
! Based on the CHIMERA code by Mike Line, but cleaned up slightly
! Pros: Fast, accurate at high optical depths, well used and familiar method
! Cons:
! IN DEVELOPMENT - USE WITH UTMOST CAUTION, YOU WILL GET WRONG ANSWERS
!!!

module ts_Toon_scatter_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Gauss quadrature variables, cosine angle values (uarr] and weights (w)
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these

  !! single angle diffusion factor approximation - typically 1/1.66
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/1.0_dp/1.66_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = uarr * w


  !! Legendre quadrature for 2 nodes
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.21132487_dp, 0.78867513_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)
  real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! Lacis & Oinas (1991) 3 point numerical values - Does not work somehow, e-mail me if you know why ;)
  ! integer, parameter :: nmu = 3
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.1_dp, 0.5_dp, 1.0_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = (/0.0433_dp, 0.5742_dp, 0.3825_dp/)

  !! Legendre quadrature for 4 nodes
  ! integer, parameter :: nmu = 4
  ! real(dp), dimension(nmu), parameter :: uarr = &
  !   & (/0.06943184_dp, 0.33000948_dp, 0.66999052_dp, 0.93056816_dp/)
  ! real(dp), dimension(nmu), parameter :: w = &
  !   & (/0.17392742_dp, 0.32607258_dp, 0.32607258_dp, 0.17392742_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = uarr * w

  !! 5 point EGP quadrature values
  ! integer, parameter :: nmu = 5
  ! real(dp), dimension(nmu), parameter :: uarr = &
  !   &(/0.0985350858_dp, 0.3045357266_dp, 0.5620251898_dp, 0.8019865821_dp, 0.9601901429_dp/)
  ! real(dp), dimension(nmu), parameter :: wuarr = &
  !   & (/0.0157479145_dp, 0.0739088701_dp, 0.1463869871_dp, 0.1671746381_dp, 0.0967815902_dp/)

  public :: ts_Toon_scatter
  private :: lw_grey_updown, sw_grey_updown, linear_log_interp

contains

  subroutine ts_Toon_scatter(nlay, nlev, Tl, pl, pe, tau_Ve, tau_IRe, mu_z, F0, Tint, AB, &
    & Beta_V, Beta_IR, sw_a, sw_g, lw_a, lw_g, net_F)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(3,nlev), intent(in) :: tau_Ve
    real(dp), dimension(2,nlev), intent(in) :: tau_IRe
    real(dp), dimension(3,nlay), intent(in) :: sw_a, sw_g
    real(dp), dimension(2,nlay), intent(in) :: lw_a, lw_g
    real(dp), dimension(3), intent(in) :: Beta_V
    real(dp), dimension(2), intent(in) :: Beta_IR
    real(dp), intent(in) :: F0, mu_z, Tint, AB

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: net_F

    !! Work variables
    integer :: i, b
    real(dp) :: be_int, Finc_b, be_int_b
    real(dp), dimension(nlev) :: Te, be, be_b
    real(dp), dimension(3,nlev) :: sw_down_b, sw_up_b
    real(dp), dimension(2,nlev) :: lw_down_b, lw_up_b
    real(dp), dimension(nlev) :: sw_down, sw_up, lw_down, lw_up
    real(dp), dimension(nlev) :: lw_net, sw_net

    !! Find temperature at layer edges through linear interpolation and extrapolation
    do i = 2, nlay
      call linear_log_interp(pe(i), pl(i-1), pl(i), Tl(i-1), Tl(i), Te(i))
      !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
    end do
    Te(1) = Tl(1) + (pe(1) - pe(2))/(pl(1) - pe(2)) * (Tl(1) - Te(2))
    Te(nlev) = Tl(nlay) + (pe(nlev) - pe(nlay))/(pl(nlay) - pe(nlay)) * (Tl(nlay) - Te(nlay))

    !! Shortwave flux calculation
    if (mu_z > 0.0_dp) then
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
      do b = 1, 3
        Finc_b = F0 * Beta_V(b)
        call sw_grey_updown(nlay, nlev, Finc_b, tau_Ve(b,:), sw_a(b,:), sw_g(b,:), mu_z, sw_up_b(b,:), sw_down_b(b,:))
        sw_down(:) = sw_down(:) + sw_down_b(b,:)
        sw_up(:) = sw_up(:) + sw_up_b(b,:)
      end do
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Longwave two-stream flux calculation
    be(:) = (sb * Te(:)**4)/pi  ! Integrated planck function intensity at levels
    be_int = (sb * Tint**4)/pi ! Integrated planck function intensity for internal temperature
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp
    do b = 1, 2
      be_b(:) = be(:) * Beta_IR(b)
      be_int_b = be_int * Beta_IR(b)
      call lw_grey_updown(nlay, nlev, be_b, be_int_b, tau_IRe(b,:), lw_up_b(b,:), lw_down_b(b,:))
      lw_up(:) = lw_up(:) + lw_up_b(b,:)
      lw_down(:) = lw_down(:) + lw_down_b(b,:)
    end do

    !! Net fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)
    sw_net(:) = sw_up(:) - sw_down(:)
    net_F(:) = lw_net(:) + sw_net(:)

    ! Uncomment if used CHIMERA lower boundary condition
    !net_F(nlev) = be_int * pi

  end subroutine ts_Toon_scatter

  subroutine lw_grey_updown(nlay, nlev, be, be_int, tau_IRe, lw_up, lw_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_IRe
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down

    !! Work variables and arrays
    integer :: k, m
    real(dp) :: tautop
    real(dp), dimension(nlay) :: dtau
    real(dp), dimension(nlay) :: B1, B0, sigma1, sigma2, em2
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    !! Calculate dtau and source functions in each layer
    do k = 1, nlay
      dtau(k) = tau_IRe(k+1) - tau_IRe(k)
      if (dtau(k) < 1e-6_dp) then
        ! For low optical depths use the isothermal approimation
        B1(k) = 0.0_dp
        B0(k) = 0.5_dp*(be(k+1) + be(k))
      else
        B1(k) = (be(k+1) - be(k))/dtau(k) ! Linear in tau term
        B0(k) = be(k)
      endif
      sigma1(k) = twopi * B0(k)
      sigma2(k) = twopi * B1(k)
    end do

    ! Upper tau boundary - (from CHIMERA)
    tautop = dtau(1)*exp(-1.0_dp)
    !tautop = 0.0_dp

    ! Zero the total flux arrays
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nmu

      !! Begin two-stream loops
      !! Perform downward loop first
      ! Top boundary condition - intensity downward from top boundary (tautop, assumed isothermal)
      lw_down_g(1) = twopi*(1.0_dp - exp(-tautop/uarr(m)))*be(1)
      do k = 1, nlay
        em2(k) = exp(-dtau(k)/uarr(m)) ! Transmission function, saved for upward loop
        lw_down_g(k+1) = lw_down_g(k)*em2(k) + sigma1(k)*(1.0_dp - em2(k)) + sigma2(k)*(uarr(m)*em2(k)+dtau(k)-uarr(m)) ! TS intensity
      end do

      !! Perform upward loop
      ! Lower boundary condition - internal heat definition Fint = F_up - F_down
      ! here we use the same condition but use intensity units to be consistent
      ! alternative boundary conditions used in CHIMERA:
      !lw_up_g(nlev) = twopi*(be(nlev) + B1(nlay)*uarr(m))
      ! if these boundary conditions are used, set net_F(nlev) = F_int (see above)
      lw_up_g(nlev) = lw_down_g(nlev) + twopi*be_int
      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*em2(k) + sigma1(k)*(1.0_dp - em2(k)) + sigma2(k)*(uarr(m)-(dtau(k)+uarr(m))*em2(k)) ! TS intensity
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      lw_down(:) = lw_down(:) + lw_down_g(:) * wuarr(m)
      lw_up(:) = lw_up(:) + lw_up_g(:) * wuarr(m)

    end do

  end subroutine lw_grey_updown

  subroutine sw_grey_updown(nlay, nlev, Finc, tau_V1, w01, gin, mu_z, sw_up, sw_down)
    implicit none

    !! Input
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: Finc, mu_z
    real(dp), dimension(nlev), intent(in) :: tau_V1
    real(dp), dimension(nlay), intent(in) :: w01, gin

    !! Output
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down

    !! Work variables
    integer :: k, i, n
    integer :: l, lm2, lm1
    real(dp) :: bsurf, rsurf, btop
    real(dp), dimension(nlev) :: direct, tau_V
    real(dp), dimension(nlay) :: dtau1, dtau
    real(dp), dimension(nlay) :: w0, g
    real(dp), dimension(nlay) :: g1, g2, g3, g4
    real(dp), dimension(nlay) :: lam, gam, alp, denom
    real(dp), dimension(nlay) :: Am, Ap, Cpm1, Cmm1, Cp, Cm
    real(dp), dimension(nlay) :: exptrm, Ep, Em, E1, E2, E3, E4
    real(dp), dimension(nlay+nlay) :: Af, Bf, Cf, Df, xk
    real(dp), dimension(nlay) :: xk1, xk2

    bsurf = 0.0_dp
    rsurf = 0.0_dp
    btop = 0.0_dp

    l = nlay + nlay
    lm2 = l - 2
    lm1 = l -1

    do k = 1, nlay
      dtau1(k) = max(tau_V1(k+1) - tau_V1(k),1e-5_dp)
    end do

    ! Delta eddington scaling
    w0(:) = (1.0_dp - gin(:)**2)*w01(:)/(1.0_dp-w01(:)*gin(:)**2)
    dtau(:) = (1.0_dp-w01(:)*gin(:)**2)*dtau1(:)
    g(:) = gin(:)/(1.0_dp + gin(:))

    ! w0(:) = w01(:)
    ! dtau(:) = dtau1(:)
    ! g(:) = gin(:)

    tau_V(1) = 0.0_dp
    do k = 1, nlay
      tau_V(k+1) = tau_V(k) + dtau(k)
    end do

    direct(:) = Finc * mu_z * exp(-tau_V(:)/mu_z)

    g1(:) = 0.86602540378_dp * (2.0_dp-w0(:)*(1.0_dp+g(:)))
    g2(:) = (1.7320508075688772_dp*w0(:)/2.0_dp) * (1.0_dp-g(:))
    g3(:) = (1.0_dp - 1.7320508075688772_dp*g(:)*mu_z)/2.0_dp
    g4(:) = 1.0_dp - g3(:)

        ! do i = 1, nlay
        !   print*, i, g1(i), g2(i), g3(i), g4(i)
        ! end do
        !
        ! stop

    lam(:) = sqrt(g1(:)**2 - g2(:)**2)
    gam(:) = (g1(:) - lam(:))/max(g2(:),1e-10_dp)
    alp(:) = sqrt((1.0_dp - w0(:))/1.0_dp - w0(:)*g(:))

    ! do i = 1, nlay
    !   print*, i, lam(i), gam(i), alp(i)
    ! end do
    !
    ! stop
    denom(:) = max(lam(:)**2 - 1.0_dp/(mu_z**2),1e-10_dp)
    Am(:) = Finc * w0(:) *(g4(:) * (g1(:) + 1.0_dp/mu_z) + g2(:)*g3(:))/denom(:)
    Ap(:) = Finc * w0(:) *(g3(:) * (g1(:) - 1.0_dp/mu_z) + g2(:)*g4(:))/denom(:)
    !
    ! do i = 1, nlay
    !   print*, i, Am(i), Ap(i)
    ! end do
    !
    ! stop


    !Cpm1 and Cmm1 are the C+ and C- terms evaluated at the top of the layer.
    Cpm1(:) = Ap(:) * exp(-tau_V(1:nlay)/mu_z)
    Cmm1(:) = Am(:) * exp(-tau_V(1:nlay)/mu_z)
    !Cp and Cm are the C+ and C- terms evaluated at the bottom of the layer.
    Cp(:) = Ap(:) * exp(-tau_V(2:nlev)/mu_z)
    Cm(:) = Am(:) * exp(-tau_V(2:nlev)/mu_z)

    ! do i = 1, nlay
    !   print*, i, Cpm1(i), Cmm1(i), Cp(i), Cm(i)
    ! end do
    !
    ! stop


    !Solve for the coefficients of system of equations using boundary conditions
    !Exponential terms:
    exptrm(:) = min(lam(:)*dtau(:),35.0_dp)
    Ep(:) = exp(exptrm(:))
    Em(:) = 1.0_dp/Ep(:)

    E1(:) = Ep(:) + gam(:)*Em(:)
    E2(:) = Ep(:) - gam(:)*Em(:)
    E3(:) = gam(:)*Ep(:) + Em(:)
    E4(:) = gam(:)*Ep(:) - Em(:)


        ! do i = 1, nlay
        !   print*, i, E1(i), E2(i), E3(i), E4(i)
        ! end do
        !
        ! stop

    Af(1) = 0.0_dp
    Bf(1) = gam(1) + 1.0_dp
    Cf(1) = gam(1) - 1.0_dp
    Df(1) = btop - Cmm1(1)

    n = 0
    do i = 2, lm2, 2
      n = n + 1
      Af(i) = (E1(n)+E3(n))*(gam(n+1)-1.0_dp)
      Bf(i) = (E2(n)+E4(n))*(gam(n+1)-1.0_dp)
      Cf(i) = 2.0_dp*(1.0_dp-gam(n+1)**2)
      Df(i) = (gam(n+1)-1.0_dp)*(Cpm1(n+1) - Cp(n)) + (1.0_dp-gam(n+1))*(Cm(n)-Cmm1(n+1))
    end do

    n = 0
    do i = 3, lm1, 2
      n = n + 1
      Af(i) = 2.0_dp*(1.0_dp-gam(n)**2)
      Bf(i) = (E1(n)-E3(n))*(1.0_dp + gam(n+1))
      Cf(i) = (E1(n)+E3(n))*(gam(n+1)-1.0_dp)
      Df(i) = E3(n)*(Cpm1(n+1) - Cp(n)) + E1(n)*(Cm(n) - Cmm1(n+1))
    end do

    Af(l) = E1(nlay) - rsurf*E3(nlay)
    Bf(l) = E2(nlay) - rsurf*E4(nlay)
    Cf(l) = 0.0_dp
    Df(l) = bsurf - Cp(nlay) + rsurf*Cm(nlay)

    ! do i = 1, l
    !   print*, i, Af(i), Bf(i), Cf(i), Df(i)
    ! end do
    !
    ! stop

    call dtridgl(l, Af, Bf, Cf, Df, xk)

    do n = 1, nlay
      xk1(n) = xk(2*n-1)+xk(2*n)
      xk2(n) = xk(2*n-1)-xk(2*n)
      if (xk2(n) == 0.0_dp) then
        cycle
      end if
      if (abs(xk2(n)/xk(2*n-1)) < 1e-30_dp) then
        xk2(n) = 0.0_dp
      end if
      !print*, xk1(n), xk2(n)
    end do

    ! sw_up(1:nlay) = xk1(:)*Ep(:)+gam(:)*xk2(:)*Em(:)+Cpm1(:)
    ! sw_down(1:nlay) = xk1(:)*Ep(:)*gam(:)+xk2(:)*Em(:)+Cmm1(:)

    sw_up(1:nlay) = xk1(:)+gam(:)*xk2(:)+Cpm1(:)
    sw_down(1:nlay) = xk1(:)*gam(:)+xk2(:)+Cmm1(:)

    sw_up(nlev) = xk1(nlay)*Ep(nlay)+gam(nlay)*xk2(nlay)*Em(nlay)+Cp(nlay)
    sw_down(nlev) = xk1(nlay)*Ep(nlay)*gam(nlay)+xk2(nlay)*Em(nlay)+Cm(nlay)

    ! do i = 1, nlev
    !   print*, i, sw_up(i), sw_down(i), direct(i)
    ! end do
    ! stop

    sw_down(:) = sw_down(:) + direct(:)
    sw_up(:) = sw_up(:)

  end subroutine sw_grey_updown


  subroutine dtridgl(l, af, bf, cf, df, xk)
    implicit none

    integer, intent(in) :: l
    real(dp), dimension(l), intent(in) :: af, bf, cf, df
    real(dp), dimension(l), intent(out) :: xk

    integer :: i
    integer, parameter :: nmax = 301
    real(dp) :: x, xkb
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

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    lxval = log10(xval)
    lx1 = log10(x1); lx2 = log10(x2)
    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / (lx2 - lx1)

    yval = 10.0_dp**((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)

  end subroutine linear_log_interp

end module ts_Toon_scatter_mod
