!!!
! Elspeth KH Lee - Jun 2024 : Initial version
!                - Oct 2026 : Rewritten from Exo-FMS_2 (fluxes at the interfaces, flux limiter, overshoot)
!
! conv: Mixing length theory (MLT) dry convection (Robinson & Marley 2014; Joyce & Tayar 2023)
!       Finds the convective heat flux at each interface where the gradient is superadiabatic,
!       and the temperature tendency from its divergence. Also returns Kzz, with convective
!       overshoot above each convective zone (Woitke et al. 2020).
!       The flux is limited so the column relaxes towards, but never past, the adiabat.
!     Pros: Time dependent convection, gives Kzz for tracer mixing, enthalpy conserving
!     Cons: Not instant (column is only near adiabatic after the convective timescale), slower than adjustment
!
! Rd and cp can be scalars or given per layer (Rd(nlay), cp(nlay)).
! Kzz and w_osh are returned at the levels (pe).
!!!

module MLT_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64
  private
  public :: MLT

  !! Rd, cp scalars or per layer. Per layer, interface k uses the mean of its two layers for
  !! Rd, cp and kappa, while each layer's heating uses its own cp (so sum cp T dp is conserved)
  interface MLT
    module procedure MLT_scalar, MLT_layer
  end interface MLT

  !! MLT parameters
  real(dp), parameter :: alpha_default = 1.0_dp   ! mixing length / pressure scale height
  real(dp), parameter :: beta_default  = 2.2_dp   ! overshoot decay exponent (Woitke et al. 2020)
  real(dp), parameter :: dt_sub_max_default = 10.0_dp  ! maximum substep [s]
  real(dp), parameter :: Kzz_min = 1.0e1_dp       ! minimum Kzz [m2/s] (1e5 cm2/s)
  real(dp), parameter :: Kzz_max = 1.0e8_dp       ! maximum Kzz [m2/s]
  real(dp), parameter :: f1 = 1.0_dp/8.0_dp, f2 = 0.5_dp
  real(dp), parameter :: limit_frac = 0.5_dp      ! fraction of the pair-neutralising flux allowed

contains

  !! Constant Rd and cp: call the per layer version with uniform arrays
  subroutine MLT_scalar(nlay, Tl, pl, pe, Rd, cp, grav, dt, dT_conv, Kzz, w_osh, alpha, beta, nsub, dt_sub_max)
    integer, intent(in) :: nlay
    real(dp), dimension(nlay),   intent(in) :: Tl, pl
    real(dp), dimension(nlay+1), intent(in) :: pe
    real(dp), intent(in) :: Rd, cp, grav, dt
    real(dp), dimension(nlay), intent(out) :: dT_conv
    real(dp), dimension(nlay+1), intent(out), optional :: Kzz
    real(dp), dimension(nlay+1), intent(out), optional :: w_osh
    real(dp), intent(in), optional :: alpha, beta, dt_sub_max
    integer, intent(in), optional :: nsub

    real(dp), dimension(nlay) :: Rd_l, cp_l

    Rd_l(:) = Rd
    cp_l(:) = cp
    call MLT_layer(nlay, Tl, pl, pe, Rd_l, cp_l, grav, dt, dT_conv, Kzz, w_osh, alpha, beta, nsub, dt_sub_max)
  end subroutine MLT_scalar

  subroutine MLT_layer(nlay, Tl, pl, pe, Rd, cp, grav, dt, dT_conv, Kzz, w_osh, alpha, beta, nsub, dt_sub_max)
    integer, intent(in) :: nlay
    real(dp), dimension(nlay),   intent(in) :: Tl, pl
    real(dp), dimension(nlay+1), intent(in) :: pe
    real(dp), dimension(nlay),   intent(in) :: Rd, cp     ! [J/kg/K] per layer
    real(dp), intent(in) :: grav, dt
    real(dp), dimension(nlay), intent(out) :: dT_conv
    real(dp), dimension(nlay+1), intent(out), optional :: Kzz     ! [m2/s] at the levels
    real(dp), dimension(nlay+1), intent(out), optional :: w_osh   ! overshoot velocity [m/s] at the levels
    real(dp), intent(in), optional :: alpha   ! mixing length / Hp (default 1)
    real(dp), intent(in), optional :: beta    ! overshoot exponent (default 2.2)
    integer, intent(in), optional :: nsub     ! minimum number of substeps (default 1)
    real(dp), intent(in), optional :: dt_sub_max   ! maximum substep [s] (default 10)

    integer :: k, isub, ns
    real(dp) :: al, be, dt_sub, Fmax, dsm
    real(dp), dimension(nlay) :: T, d_p, dT_s, kap_l
    real(dp), dimension(nlay+1) :: Fc, F_mlt, w, grad, T_e, Kzz_e, w_o, Kzz_acc, w_acc
    real(dp), dimension(nlay+1) :: Rd_e, cp_e, kap_e   ! at the interfaces 2..nlay
    logical :: want_kzz
    real(dp), dimension(nlay-1) :: pfact
    logical, dimension(nlay+1) :: conv

    al = alpha_default;  if (present(alpha)) al = alpha
    be = beta_default;   if (present(beta)) be = beta
    dsm = dt_sub_max_default;  if (present(dt_sub_max)) dsm = dt_sub_max

    !! Number of substeps (at most dt_sub_max each, at least nsub)
    ns = max(1, ceiling(dt/dsm))
    if (present(nsub)) ns = max(ns, nsub)
    dt_sub = dt/real(ns, dp)

    !! Gas properties at the interfaces - mean of the two layers
    kap_l(:) = Rd(:)/cp(:)
    Rd_e(:) = Rd(1);  cp_e(:) = cp(1);  kap_e(:) = kap_l(1)
    Rd_e(2:nlay)  = 0.5_dp*(Rd(1:nlay-1) + Rd(2:nlay))
    cp_e(2:nlay)  = 0.5_dp*(cp(1:nlay-1) + cp(2:nlay))
    kap_e(2:nlay) = 0.5_dp*(kap_l(1:nlay-1) + kap_l(2:nlay))

    !! Layer pressure thickness and neutral pair factor: T(k-1) = T(k) pfact(k-1) is neutral
    d_p(:) = pe(2:nlay+1) - pe(1:nlay)
    pfact(:) = (pl(1:nlay-1)/pl(2:nlay))**kap_e(2:nlay)
    T(:) = Tl(:)
    want_kzz = present(Kzz) .or. present(w_osh)
    Kzz_acc(:) = 0.0_dp;  w_acc(:) = 0.0_dp

    !! Explicit substeps
    do isub = 1, ns
      call interface_state(T, grad, T_e, conv, w, F_mlt)

      !! Limit each flux to limit_frac of the flux that would mix its pair to neutral in dt_sub
      Fc(:) = F_mlt(:)
      do k = 2, nlay
        if (.not. conv(k)) cycle
        Fmax = (cp(k-1)/grav)*d_p(k-1)*(pair_neutral_upper(k) - T(k-1))/dt_sub
        Fc(k) = min(Fc(k), limit_frac*max(Fmax, 0.0_dp))
      end do

      !! Kzz and overshoot of this substep, with the velocity of the flux actually carried (F ~ w^3)
      if (want_kzz) then
        do k = 2, nlay
          if (conv(k) .and. F_mlt(k) > 0.0_dp) w(k) = w(k)*(Fc(k)/F_mlt(k))**(1.0_dp/3.0_dp)
        end do
        call kzz_overshoot(conv, w, T_e, Kzz_e, w_o)
        Kzz_acc(:) = Kzz_acc(:) + Kzz_e(:)
        w_acc(:) = w_acc(:) + w_o(:)
      end if

      !! Heating from the flux divergence
      do k = 1, nlay
        dT_s(k) = (grav/cp(k))*(Fc(k+1) - Fc(k))/d_p(k)
      end do
      T(:) = T(:) + dt_sub*dT_s(:)
    end do

    dT_conv(:) = (T(:) - Tl(:))/dt

    !! Kzz and overshoot velocity averaged over the substeps
    if (present(Kzz)) Kzz(:) = min(Kzz_acc(:)/real(ns, dp), Kzz_max)
    if (present(w_osh)) w_osh(:) = w_acc(:)/real(ns, dp)

  contains

    !! Kzz at the interfaces (w L where convective, Kzz_min elsewhere) and overshoot (L w_osh)
    !! above the top of each convective zone, up to the next zone above
    subroutine kzz_overshoot(conv_i, w_i, T_i, Kzz_i, wo_i)
      logical, dimension(nlay+1), intent(in) :: conv_i
      real(dp), dimension(nlay+1), intent(in)  :: w_i, T_i
      real(dp), dimension(nlay+1), intent(out) :: Kzz_i, wo_i
      integer :: i, it
      real(dp) :: Hp, wov
      Kzz_i(:) = Kzz_min
      wo_i(:) = 0.0_dp
      do i = 2, nlay
        if (conv_i(i)) Kzz_i(i) = max(w_i(i)*al*Rd_e(i)*T_i(i)/grav, Kzz_min)
      end do
      do it = 2, nlay
        if (.not. conv_i(it)) cycle
        if (conv_i(it-1)) cycle                            ! not the top of its zone
        do i = it-1, 2, -1
          if (conv_i(i)) exit
          Hp = Rd_e(i)*T_i(i)/grav
          wov = w_i(it)*(pe(i)/pe(it))**be
          wo_i(i) = max(wo_i(i), wov)
          Kzz_i(i) = max(Kzz_i(i), al*Hp*wov)
        end do
      end do
      ! Bottom level has no interface below it, use the value of the level above
      Kzz_i(nlay+1) = Kzz_i(nlay)
    end subroutine kzz_overshoot

    !! Gradient, temperature, convective flag, velocity and MLT flux at the interfaces
    subroutine interface_state(Tc, grad_i, T_i, conv_i, w_i, F_i)
      real(dp), dimension(nlay),   intent(in)  :: Tc
      real(dp), dimension(nlay+1), intent(out) :: grad_i, T_i, w_i, F_i
      logical, dimension(nlay+1), intent(out) :: conv_i
      integer :: i
      real(dp) :: Hp, L, rho
      grad_i(:) = 0.0_dp;  w_i(:) = 0.0_dp;  F_i(:) = 0.0_dp;  conv_i(:) = .false.
      T_i(1) = Tc(1);  T_i(nlay+1) = Tc(nlay)
      do i = 2, nlay
        grad_i(i) = log(Tc(i-1)/Tc(i))/log(pl(i-1)/pl(i))
        T_i(i) = Tc(i-1)*(pe(i)/pl(i-1))**grad_i(i)      ! power law through both layers
        if (grad_i(i) > kap_e(i)) then
          conv_i(i) = .true.
          Hp = Rd_e(i)*T_i(i)/grav
          L = al*Hp
          rho = pe(i)/(Rd_e(i)*T_i(i))
          w_i(i) = sqrt(f1*L**2*grav*(grad_i(i) - kap_e(i))/Hp)
          F_i(i) = f2*rho*cp_e(i)*T_i(i)*w_i(i)*al*(grad_i(i) - kap_e(i))
        end if
      end do
    end subroutine interface_state

    !! Temperature of layer i-1 after mixing layers i-1 and i to neutral, conserving the
    !! enthalpy of the pair (r = cp(i)/cp(i-1), 1 for constant cp)
    real(dp) function pair_neutral_upper(i)
      integer, intent(in) :: i
      real(dp) :: Tlow, r
      r = cp(i)/cp(i-1)
      Tlow = (d_p(i-1)*T(i-1) + r*d_p(i)*T(i))/(r*d_p(i) + pfact(i-1)*d_p(i-1))
      pair_neutral_upper = Tlow*pfact(i-1)
    end function pair_neutral_upper

  end subroutine MLT_layer

end module MLT_mod
