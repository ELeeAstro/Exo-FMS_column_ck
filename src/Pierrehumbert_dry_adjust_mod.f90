!!!
! Elspeth KH Lee - Sep 2020 : Initial version
!                - Oct 2026 : Updated from Exo-FMS_2 (tracer and wind mixing, convergence flag)
!
! conv: Dry convective adjustment following Ray Pierrehumbert's python notebook implementation
!       Sweeps down and then up the column, mixing each unstable pair of layers to neutral while
!       conserving the enthalpy of the pair, repeated until no pair needs adjusting.
!     Pros: Simple and robust, same profile as the Manabe scheme when converged
!     Cons: A correction moves one layer per pass, so deep unstable regions need many passes
!           (the column can be left unstable if itermax is reached, see converged)
!
! Optional: tracers (q) and winds (u, v) are mixed to their mass-weighted pair mean as each pair
!           is adjusted, with the kinetic energy lost by the wind mixing returned as heat if cp is given.
! Pass Tl including the other tendencies of the step (T + dt dT_other), so the column is
! stable at the end of the step.
!!!

module Pierrehumbert_dry_adjust_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64
  private
  public :: Pierrehumbert_dry_adjust

  !! Maximum number of down + up passes and instability threshold [K]
  integer, parameter :: itermax = 100
  real(dp),    parameter :: small   = 1.0e-6_dp

contains

  subroutine Pierrehumbert_dry_adjust(nlay, Tl, pl, pe, kappa, dt, dT_conv, niter, converged, &
    & q, dq_conv, u, v, du_conv, dv_conv, cp)
    integer, intent(in) :: nlay
    real(dp), dimension(nlay),   intent(in) :: Tl, pl
    real(dp), dimension(nlay+1), intent(in) :: pe
    real(dp), intent(in) :: kappa, dt
    real(dp), dimension(nlay), intent(out) :: dT_conv
    integer, intent(out), optional :: niter      ! number of passes used
    logical, intent(out), optional :: converged  ! .false. if itermax was reached
    real(dp), dimension(:,:), intent(in),  optional :: q        ! tracers (nlay, nq)
    real(dp), dimension(:,:), intent(out), optional :: dq_conv  ! tracer tendencies (nlay, nq) [1/s]
    real(dp), dimension(nlay), intent(in),  optional :: u, v    ! winds [m/s]
    real(dp), dimension(nlay), intent(out), optional :: du_conv, dv_conv  ! [m/s2]
    real(dp), intent(in), optional :: cp          ! [J/kg/K] return the wind-mixing KE loss as heat

    integer :: i, iter
    logical :: did_adj, mix_q, mix_uv
    real(dp), dimension(nlay) :: Tl_cc, d_p, uc, vc
    real(dp), dimension(nlay-1) :: pfact
    real(dp), allocatable :: qc(:,:)
    real(dp) :: Tbar

    mix_q = present(q) .and. present(dq_conv)
    mix_uv = present(u) .and. present(v) .and. present(du_conv) .and. present(dv_conv)
    if (present(q) .neqv. present(dq_conv)) error stop 'Pierrehumbert_dry_adjust: give q and dq_conv together'
    if ((present(u) .or. present(v) .or. present(du_conv) .or. present(dv_conv)) .and. .not. mix_uv) &
      & error stop 'Pierrehumbert_dry_adjust: give u, v, du_conv and dv_conv together'

    !! Work on a copy of the temperature (and tracers and winds)
    Tl_cc(:) = Tl(:)
    if (mix_q) qc = q
    if (mix_uv) then
      uc(:) = u(:);  vc(:) = v(:)
    end if

    !! Layer pressure thickness and neutral pair factor: T(i) = T(i+1) pfact(i) is neutral
    d_p(:) = pe(2:nlay+1) - pe(1:nlay)
    pfact(:) = (pl(1:nlay-1)/pl(2:nlay))**kappa

    did_adj = .true.
    do iter = 1, itermax

      did_adj = .false.

      !! Downward pass
      do i = 1, nlay-1
        if (Tl_cc(i) < (Tl_cc(i+1)*pfact(i) - small)) then
          call adjust_pair(i)
        end if
      end do

      !! Upward pass
      do i = nlay-1, 1, -1
        if (Tl_cc(i) < (Tl_cc(i+1)*pfact(i) - small)) then
          call adjust_pair(i)
        end if
      end do

      !! If no adjustment required, exit the loop
      if (.not. did_adj) exit

    end do

    if (present(niter)) niter = min(iter, itermax)
    if (present(converged)) converged = .not. did_adj

    !! Change in temperature is Tl_cc - Tl
    !! adjust on timescale of 1 timestep (i.e. instant correction across one timestep)
    dT_conv(:) = (Tl_cc(:) - Tl(:))/dt
    if (mix_q) dq_conv(:,:) = (qc(:,:) - q(:,:))/dt
    if (mix_uv) then
      du_conv(:) = (uc(:) - u(:))/dt;  dv_conv(:) = (vc(:) - v(:))/dt
    end if

  contains

    !! Mix layers i and i+1 to neutral, conserving the enthalpy of the pair,
    !! and tracers and winds to their mass-weighted pair mean
    subroutine adjust_pair(i)
      integer, intent(in) :: i
      real(dp) :: m1, m2, um, vm, dke, c
      m1 = d_p(i);  m2 = d_p(i+1)
      Tbar = (m1*Tl_cc(i) + m2*Tl_cc(i+1))/(m1 + m2)
      Tl_cc(i+1) = (m1 + m2)*Tbar/(m2 + pfact(i)*m1)
      Tl_cc(i) = Tl_cc(i+1)*pfact(i)
      if (mix_q) then
        qc(i,:) = (m1*qc(i,:) + m2*qc(i+1,:))/(m1 + m2)
        qc(i+1,:) = qc(i,:)
      end if
      if (mix_uv) then
        um = (m1*uc(i) + m2*uc(i+1))/(m1 + m2)
        vm = (m1*vc(i) + m2*vc(i+1))/(m1 + m2)
        dke = 0.5_dp*(m1*(uc(i)**2 + vc(i)**2) + m2*(uc(i+1)**2 + vc(i+1)**2) - (m1 + m2)*(um**2 + vm**2))
        uc(i) = um;  uc(i+1) = um;  vc(i) = vm;  vc(i+1) = vm
        if (present(cp)) then
          ! KE lost -> heat, keeping the pair neutral: dT(i+1) = c, dT(i) = c pfact(i)
          c = dke/(cp*(m2 + pfact(i)*m1))
          Tl_cc(i+1) = Tl_cc(i+1) + c
          Tl_cc(i) = Tl_cc(i) + c*pfact(i)
        end if
      end if
      did_adj = .true.
    end subroutine adjust_pair

  end subroutine Pierrehumbert_dry_adjust

end module Pierrehumbert_dry_adjust_mod
