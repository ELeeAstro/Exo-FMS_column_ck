!!!
! Elspeth KH Lee - Oct 2026 : Initial version (from Exo-FMS_2)
!
! conv: Dry convective adjustment following Manabe, Smagorinsky & Strickler (1965)
!       Each unstable segment of the column is set to a single potential temperature, conserving
!       the enthalpy of the segment. A single bottom-up sweep merges unstable layers into segments.
!       Gives the same profile as the Pierrehumbert scheme when that one is fully converged.
!     Pros: No iteration (cost linear in nlay), exactly neutral and enthalpy conserving
!     Cons: Instant adjustment over one timestep (use relax_frac for partial relaxation)
!
! Optional: tracers (q) and winds (u, v) are mixed to their mass-weighted mean in each segment,
!           with the kinetic energy lost by the wind mixing returned as heat if cp is given.
! Pass Tl including the other tendencies of the step (T + dt dT_other), so the column is
! stable at the end of the step.
!!!

module Manabe_dry_adjust_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64
  private
  public :: Manabe_dry_adjust

contains

  subroutine Manabe_dry_adjust(nlay, Tl, pl, pe, kappa, dt, dT_conv, nadj, &
    & q, dq_conv, u, v, du_conv, dv_conv, cp, relax_frac)
    integer, intent(in) :: nlay
    real(dp), dimension(nlay),   intent(in) :: Tl, pl
    real(dp), dimension(nlay+1), intent(in) :: pe
    real(dp), intent(in) :: kappa, dt
    real(dp), dimension(nlay), intent(out) :: dT_conv
    integer, intent(out), optional :: nadj   ! number of layers in adjusted segments
    real(dp), dimension(:,:), intent(in),  optional :: q        ! tracers (nlay, nq)
    real(dp), dimension(:,:), intent(out), optional :: dq_conv  ! tracer tendencies (nlay, nq) [1/s]
    real(dp), dimension(nlay), intent(in),  optional :: u, v    ! winds [m/s]
    real(dp), dimension(nlay), intent(out), optional :: du_conv, dv_conv  ! [m/s2]
    real(dp), intent(in), optional :: cp          ! [J/kg/K] return the wind-mixing KE loss as heat
    real(dp), intent(in), optional :: relax_frac  ! fraction of the adjustment applied (default 1)

    integer :: k, m, n, nseg
    integer, dimension(nlay) :: seg_top, seg_bot
    real(dp), dimension(nlay) :: seg_T, seg_P    ! sum T dp and sum Pi dp of each segment
    real(dp), dimension(nlay) :: Pi, d_p
    integer :: top, bot
    real(dp) :: sT, sP, f, mass, qm, um, vm, dke, unew, vnew
    logical :: mix_q, mix_uv

    mix_q = present(q) .and. present(dq_conv)
    mix_uv = present(u) .and. present(v) .and. present(du_conv) .and. present(dv_conv)
    if (present(q) .neqv. present(dq_conv)) error stop 'Manabe_dry_adjust: give q and dq_conv together'
    if ((present(u) .or. present(v) .or. present(du_conv) .or. present(dv_conv)) .and. .not. mix_uv) &
      & error stop 'Manabe_dry_adjust: give u, v, du_conv and dv_conv together'
    f = 1.0_dp
    if (present(relax_frac)) f = relax_frac

    !! Layer pressure thickness and Exner function (the reference pressure cancels)
    d_p(:) = pe(2:nlay+1) - pe(1:nlay)
    Pi(:) = (pl(:)/pl(nlay))**kappa

    !! Bottom-up sweep keeping a stack of neutral segments (seg_*(nseg) is the one just below)
    nseg = 0
    do k = nlay, 1, -1
      top = k;  bot = k
      sT = Tl(k)*d_p(k);  sP = Pi(k)*d_p(k)
      do while (nseg > 0)
        if (sT*seg_P(nseg) >= seg_T(nseg)*sP) exit    ! theta(new) >= theta(below): stable
        bot = seg_bot(nseg)                             ! unstable: merge with the segment below
        sT = sT + seg_T(nseg);  sP = sP + seg_P(nseg)
        nseg = nseg - 1
      end do
      nseg = nseg + 1
      seg_top(nseg) = top;  seg_bot(nseg) = bot
      seg_T(nseg) = sT;     seg_P(nseg) = sP
    end do

    !! Neutral profile (and mixed tracers and winds) in each merged segment, single layers unchanged
    dT_conv(:) = 0._dp
    if (mix_q) dq_conv(:,:) = 0._dp
    if (mix_uv) then
      du_conv(:) = 0._dp;  dv_conv(:) = 0._dp
    end if
    if (present(nadj)) nadj = 0
    do n = 1, nseg
      top = seg_top(n);  bot = seg_bot(n)
      if (bot == top) cycle
      do k = top, bot
        dT_conv(k) = f*(seg_T(n)/seg_P(n)*Pi(k) - Tl(k))/dt
      end do
      mass = sum(d_p(top:bot))
      if (mix_q) then
        do m = 1, size(q, 2)
          qm = sum(q(top:bot,m)*d_p(top:bot))/mass
          dq_conv(top:bot,m) = f*(qm - q(top:bot,m))/dt
        end do
      end if
      if (mix_uv) then
        um = sum(u(top:bot)*d_p(top:bot))/mass
        vm = sum(v(top:bot)*d_p(top:bot))/mass
        dke = 0._dp
        do k = top, bot
          unew = u(k) + f*(um - u(k));  vnew = v(k) + f*(vm - v(k))
          du_conv(k) = (unew - u(k))/dt;  dv_conv(k) = (vnew - v(k))/dt
          dke = dke + 0.5_dp*(u(k)**2 + v(k)**2 - unew**2 - vnew**2)*d_p(k)   ! >= 0 summed
        end do
        if (present(cp)) then
          ! KE lost -> heat, uniform in potential temperature over the segment
          dT_conv(top:bot) = dT_conv(top:bot) + dke/(cp*seg_P(n))*Pi(top:bot)/dt
        end if
      end if
      if (present(nadj)) nadj = nadj + bot - top + 1
    end do

  end subroutine Manabe_dry_adjust

end module Manabe_dry_adjust_mod
