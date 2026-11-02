!!!
! Elspeth KH Lee - Nov 2026 : Initial version, converted from AddingDoublingRT (JAX/C++)
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!
! sw: Adding-doubling method (matrix operator method; Plass, Kattawar & Catchings 1973;
!     Hansen 1971), converted to Fortran from AddingDoublingRT (src_jax/doubling.py,
!     adding.py, solver.py; same algorithm as the C++ src/), stellar beam source only.
!     Doubling:  reflection R, transmission T and beam source vectors of each homogeneous
!                layer, started from an optically thin sub-layer (single scattering of the
!                beam) and doubled nn times.
!     Adding:    one-sided (Riccati): composites from the bottom (reflection for light from above, RbR,
!                and upward beam source, Rbs - the only parts needed), then a single downward sweep
!                I_dn(l+1) = (I - R_l RbR_{l+1})^-1 [T_l I_dn(l) + sdn_l + R_l Rbs_{l+1}],
!                I_up(l) = Rbs_l + RbR_l I_dn(l). Identical to the two-sided composites (rtop/rbase).
!     Follows the zcorr beam path (each layer has its own zenith angle) and conserves energy.
!
! Pros: accurate multiple scattering for any w and g, converges with nmu
! Cons: nn ~ log2(dtau) + 16 doublings (N x N matrix solves) per layer - the most expensive
!       sw method apart from sw_disort
!!!

module sw_AD_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

  !! Number of Gauss-Legendre streams per hemisphere (2*nmu stream method)
  integer, parameter :: nmu = 4

  !! Number of phase function Legendre moments retained (delta-M truncates chi_(2 nmu))
  integer, parameter :: nmom = 2*nmu

  !! Layer (or composite) reflection, transmission and beam sources.
  !! ab: incident from the top (a) side, ba: incident from the bottom (b) side.
  type :: ad_layer
    real(dp), dimension(nmu,nmu) :: Rab, Rba, Tab, Tba
    real(dp), dimension(nmu) :: sup, sdn
  end type ad_layer

  private :: sw_adding_doubling, layer_doubling, gauss_legendre, &
    & right_div, left_solve, lu_factor, lu_solve, ipow0
  public :: sw_AD

contains

  subroutine sw_AD(nlay, nlev, nb, ng, gw, tau_e, mu_z, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev, nb, ng
    real(dp), dimension(ng), intent(in) :: gw
    real(dp), dimension(ng,nb,nlev), intent(in) :: tau_e
    real(dp), dimension(ng,nb,nlay), intent(in) :: ssa, gg
    real(dp), dimension(nlev), intent(in) :: mu_z
    real(dp), dimension(nb), intent(in) :: Finc

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down, sw_net
    real(dp), intent(out) :: asr

    !! Work variables
    integer :: b, g
    real(dp), dimension(nlev) :: sw_up_g, sw_down_g

    sw_up(:) = 0.0_dp
    sw_down(:) = 0.0_dp

    !! Shortwave flux calculation - sum over each band and g-ordinance
    if (mu_z(nlev) > 0.0_dp) then
      do b = 1, nb
        if (Finc(b) <= 0.0_dp) cycle
        do g = 1, ng
          call sw_adding_doubling(nlay, nlev, Finc(b), mu_z(:), tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), &
            & sw_down_g(:), sw_up_g(:))
          sw_up(:) = sw_up(:) + sw_up_g(:) * gw(g)
          sw_down(:) = sw_down(:) + sw_down_g(:) * gw(g)
        end do
      end do
    end if

    !! Net sw flux
    sw_net(:) = sw_up(:) - sw_down(:)

    !! Absorbed Stellar Radiation (ASR)
    asr = sw_down(1) - sw_up(1)

  end subroutine sw_AD

  subroutine sw_adding_doubling(nlay, nlev, F0_in, mu_in, tau_in, w_in, g_in, flx_down, flx_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0_in
    real(dp), dimension(nlev), intent(in) :: tau_in, mu_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    !! Work variables and arrays
    integer :: k, l, i, j
    real(dp), dimension(nmu) :: mu, wt, Iup, Idn, ppl, pmn
    real(dp), dimension(nmom, nmu) :: Pl
    real(dp), dimension(nmom) :: chi, c2l, sgn, Pmu0
    real(dp), dimension(nmu,nmu) :: Ppp, Ppm, Imat
    real(dp), dimension(nlay) :: dtau, w0, mu_l, F_top
    real(dp), dimension(nlev) :: cum_trans, T_beam
    real(dp) :: f, colsum
    type(ad_layer), dimension(nlay) :: lay
    real(dp), dimension(nmu,nmu,nlev) :: RbR
    real(dp), dimension(nmu,nlev) :: Rbs
    real(dp), dimension(nmu,nmu) :: TD1

    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w_in(:) <= 1.0e-12_dp)) then

      ! Check mu_z for spherical corrections
      if (mu_in(nlev) == mu_in(1)) then
        ! No zenith correction, use regular method
        flx_down(:) = F0_in * mu_in(nlev) * exp(-tau_in(:)/mu_in(nlev))
      else
        ! Zenith angle correction, use cumulative transmission function
        cum_trans(1) = tau_in(1)/mu_in(1)
        do k = 1, nlev-1
          cum_trans(k+1) = cum_trans(k) + (tau_in(k+1) - tau_in(k))/mu_in(k+1)
        end do
        flx_down(:) = F0_in * mu_in(nlev) * exp(-cum_trans(:))
      end if

      flx_up(:) = 0.0_dp

      return

    end if

    !! Gauss-Legendre nodes and weights on [0,1] (sum(wt) = 1, sum(wt*mu) = 1/2)
    call gauss_legendre(nmu, mu, wt)

    !! Legendre polynomials P_l(mu_i), l = 0..nmom-1 (stored at l+1)
    Pl(1,:) = 1.0_dp
    Pl(2,:) = mu(:)
    do l = 2, nmom-1
      Pl(l+1,:) = (real(2*l-1,dp)*mu(:)*Pl(l,:) - real(l-1,dp)*Pl(l-1,:))/real(l,dp)
    end do

    !! (2l+1) and (-1)^l factors of the phase matrix expansion
    do l = 0, nmom-1
      c2l(l+1) = real(2*l+1,dp)
      sgn(l+1) = (-1.0_dp)**l
    end do

    Imat(:,:) = 0.0_dp
    do i = 1, nmu
      Imat(i,i) = 1.0_dp
    end do

    !! Layer optical depths and delta-M scaling (HG, f = chi_(2 nmu) = g**(2 nmu))
    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)
    w0(:) = min(max(w_in(:), 0.0_dp), 1.0_dp)
    do k = 1, nlay
      f = g_in(k)**nmom
      if (w0(k) > 0.0_dp .and. f > 1.0e-12_dp .and. f < 1.0_dp - 1.0e-12_dp) then
        dtau(k) = (1.0_dp - w0(k)*f)*dtau(k)
        w0(k) = w0(k)*(1.0_dp - f)/(1.0_dp - w0(k)*f)
      end if
    end do

    !! Zenith angle in each layer (zcorr uses the lower level, as the cumulative transmission)
    if (mu_in(nlev) == mu_in(1)) then
      mu_l(:) = mu_in(nlev)
    else
      mu_l(:) = mu_in(2:nlev)
    end if

    !! Direct beam transmission to each level along the beam path (scaled tau)
    cum_trans(1) = 0.0_dp
    do k = 1, nlay
      cum_trans(k+1) = cum_trans(k) + dtau(k)/mu_l(k)
    end do
    T_beam(:) = exp(-min(cum_trans(:), 200.0_dp))

    !! Beam flux (normal to the beam) at the top of each layer, so the beam scatters
    !! exactly the energy it loses (energy conserving with zcorr)
    F_top(:) = F0_in * mu_in(nlev)/mu_l(:) * T_beam(1:nlay)

    !! Doubling - R, T and beam sources of each layer
    do k = 1, nlay

      Ppp(:,:) = 0.0_dp
      Ppm(:,:) = 0.0_dp
      ppl(:) = 0.0_dp
      pmn(:) = 0.0_dp

      if (w0(k) > 0.0_dp .and. dtau(k) > 0.0_dp) then

        !! Henyey-Greenstein moments, delta-M scaled with the first truncated moment chi_(2 nmu)
        f = g_in(k)**nmom
        do l = 0, nmom-1
          chi(l+1) = g_in(k)**l
        end do
        if (f > 1.0e-12_dp .and. f < 1.0_dp - 1.0e-12_dp) then
          chi(:) = (chi(:) - f)/(1.0_dp - f)
        end if

        !! Azimuthally averaged phase matrices for same (Ppp) and opposite (Ppm) hemispheres:
        !! P(+/-mu_i, mu_j) = sum_l (2l+1) chi_l P_l(mu_i) P_l(mu_j) (+/-1)^l / (2 pi)
        do j = 1, nmu
          do i = 1, nmu
            Ppp(i,j) = sum(c2l(:)*chi(:)*Pl(:,i)*Pl(:,j))/(2.0_dp*pi)
            Ppm(i,j) = sum(sgn(:)*c2l(:)*chi(:)*Pl(:,i)*Pl(:,j))/(2.0_dp*pi)
          end do
        end do

        !! Hansen normalisation: each incident stream j scatters exactly w0 of its energy
        do j = 1, nmu
          colsum = sum((Ppp(:,j) + Ppm(:,j))*wt(:))
          if (colsum > 0.0_dp) then
            Ppp(:,j) = Ppp(:,j)/(2.0_dp*pi*colsum)
            Ppm(:,j) = Ppm(:,j)/(2.0_dp*pi*colsum)
          end if
        end do

        !! Beam phase vectors: scattering from the beam (mu0 = mu_l) into the downward (ppl)
        !! and upward (pmn) streams, normalised in the same way
        Pmu0(1) = 1.0_dp
        Pmu0(2) = mu_l(k)
        do l = 2, nmom-1
          Pmu0(l+1) = (real(2*l-1,dp)*mu_l(k)*Pmu0(l) - real(l-1,dp)*Pmu0(l-1))/real(l,dp)
        end do
        do i = 1, nmu
          ppl(i) = sum(c2l(:)*chi(:)*Pl(:,i)*Pmu0(:))/(2.0_dp*pi)
          pmn(i) = sum(sgn(:)*c2l(:)*chi(:)*Pl(:,i)*Pmu0(:))/(2.0_dp*pi)
        end do
        colsum = sum((ppl(:) + pmn(:))*wt(:))
        if (colsum > 0.0_dp) then
          ppl(:) = ppl(:)/(2.0_dp*pi*colsum)
          pmn(:) = pmn(:)/(2.0_dp*pi*colsum)
        end if

      end if

      call layer_doubling(dtau(k), w0(k), mu_l(k), F_top(k), Ppp, Ppm, ppl, pmn, mu, wt, Imat, lay(k))

    end do

    !! Upward pass - composite of everything below each level: reflection for light from above (RbR)
    !! and upward beam source (Rbs). Lower boundary: gas giant, black - nothing is reflected back up
    RbR(:,:,nlev) = 0.0_dp
    Rbs(:,nlev) = 0.0_dp
    do l = nlay, 1, -1
      TD1(:,:) = right_div(lay(l)%Tba, Imat - matmul(RbR(:,:,l+1), lay(l)%Rba))
      RbR(:,:,l) = lay(l)%Rab(:,:) + matmul(matmul(TD1, RbR(:,:,l+1)), lay(l)%Tab)
      Rbs(:,l) = lay(l)%sup(:) + matmul(TD1, Rbs(:,l+1) + matmul(RbR(:,:,l+1), lay(l)%sdn))
    end do

    !! Top of atmosphere - no diffuse intensity from above, only the beam
    Idn(:) = 0.0_dp
    flx_up(1) = 2.0_dp*pi*sum(wt(:)*mu(:)*Rbs(:,1))
    flx_down(1) = F0_in*mu_in(nlev)*T_beam(1)

    !! Downward pass: I_dn(l+1) = (I - Rba_l RbR_{l+1})^-1 [Tab_l I_dn(l) + sdn_l + Rba_l Rbs_{l+1}]
    do l = 1, nlay
      Idn(:) = left_solve(Imat - matmul(lay(l)%Rba, RbR(:,:,l+1)), &
        & matmul(lay(l)%Tab, Idn) + lay(l)%sdn + matmul(lay(l)%Rba, Rbs(:,l+1)))
      Iup(:) = Rbs(:,l+1) + matmul(RbR(:,:,l+1), Idn)
      flx_up(l+1) = 2.0_dp*pi*sum(wt(:)*mu(:)*Iup(:))
      flx_down(l+1) = 2.0_dp*pi*sum(wt(:)*mu(:)*Idn(:)) + F0_in*mu_in(nlev)*T_beam(l+1)
    end do

  end subroutine sw_adding_doubling

  !! R, T and beam sources of one homogeneous layer by doubling
  subroutine layer_doubling(tau, omega, mu0, F_top, Ppp, Ppm, ppl, pmn, mu, wt, Imat, lyr)
    implicit none

    real(dp), intent(in) :: tau, omega, mu0, F_top
    real(dp), dimension(nmu,nmu), intent(in) :: Ppp, Ppm, Imat
    real(dp), dimension(nmu), intent(in) :: ppl, pmn, mu, wt
    type(ad_layer), intent(out) :: lyr

    integer :: i, n, nn
    real(dp) :: tau0, con, gam
    real(dp), dimension(nmu) :: trans, sup, sdn, rhs_up, rhs_dn
    real(dp), dimension(nmu,nmu) :: R, T, TG, TGR, Gpp, Gpm

    lyr%Rab(:,:) = 0.0_dp
    lyr%Rba(:,:) = 0.0_dp
    lyr%Tab(:,:) = Imat(:,:)
    lyr%Tba(:,:) = Imat(:,:)
    lyr%sup(:) = 0.0_dp
    lyr%sdn(:) = 0.0_dp

    !! Zero optical depth - transparent layer
    if (tau <= 0.0_dp) then
      return
    end if

    !! Pure absorption - exact diagonal transmission, no scattered beam
    if (omega <= 0.0_dp) then
      trans(:) = exp(max(-tau/mu(:), -200.0_dp))
      lyr%Tab(:,:) = 0.0_dp
      do i = 1, nmu
        lyr%Tab(i,i) = trans(i)
      end do
      lyr%Tba(:,:) = lyr%Tab(:,:)
      return
    end if

    !! Scattering layer - single-scattering operators of an infinitesimal layer:
    !! Gpp = (I - 2 pi w Ppp C)/mu, Gpm = 2 pi w Ppm C/mu, C = diag(wt)
    con = 2.0_dp*pi*omega
    do i = 1, nmu
      Gpp(i,:) = (Imat(i,:) - con*Ppp(i,:)*wt(:))/mu(i)
      Gpm(i,:) = con*Ppm(i,:)*wt(:)/mu(i)
    end do

    !! Number of doublings: start from tau0 = tau/2^nn, thin enough for single scattering
    nn = max(int(log(tau)/log(2.0_dp)) + ipow0(omega), 1)
    tau0 = tau/(2.0_dp**nn)

    !! Initial thin sub-layer: R, T and single scattering of the beam into the streams
    R(:,:) = tau0*Gpm(:,:)
    T(:,:) = Imat(:,:) - tau0*Gpp(:,:)
    sup(:) = omega*tau0/mu(:)*F_top*pmn(:)
    sdn(:) = omega*tau0/mu(:)*F_top*ppl(:)

    !! Beam transmission across the sub-layer
    gam = exp(-tau0/mu0)

    !! Doubling: combine two identical sub-layers into one of twice the depth, the beam
    !! reaching the lower sub-layer is reduced by gam
    do n = 1, nn
      TG(:,:) = right_div(T, Imat - matmul(R, R))  ! T (I - R^2)^-1
      TGR(:,:) = matmul(TG, R)

      rhs_up(:) = matmul(R, sdn) + gam*sup(:)
      rhs_dn(:) = gam*matmul(R, sup) + sdn(:)
      sup(:) = matmul(TG, rhs_up) + sup(:)
      sdn(:) = matmul(TG, rhs_dn) + gam*sdn(:)
      gam = gam*gam

      R(:,:) = R(:,:) + matmul(TGR, T)
      T(:,:) = matmul(TG, T)
    end do

    !! Homogeneous layer: symmetric about its mid-plane
    lyr%Rab(:,:) = R(:,:)
    lyr%Rba(:,:) = R(:,:)
    lyr%Tab(:,:) = T(:,:)
    lyr%Tba(:,:) = T(:,:)
    lyr%sup(:) = sup(:)
    lyr%sdn(:) = sdn(:)

  end subroutine layer_doubling

  !! Adaptive number of initial doublings: single-scattering error is O(tau0^2 w^2)
  integer function ipow0(omega)
    real(dp), intent(in) :: omega
    if (omega < 0.01_dp) then
      ipow0 = 4
    else if (omega < 0.1_dp) then
      ipow0 = 10
    else
      ipow0 = 16
    end if
  end function ipow0

  !! Gauss-Legendre nodes and weights on [0,1]
  subroutine gauss_legendre(n, x, wq)
    implicit none

    integer, intent(in) :: n
    real(dp), dimension(n), intent(out) :: x, wq

    integer :: i, k, it
    real(dp) :: z, dz, p0, p1, pk, pp

    do i = 1, (n + 1)/2
      z = cos(pi*(real(i,dp) - 0.25_dp)/(real(n,dp) + 0.5_dp))
      do it = 1, 100
        p0 = 1.0_dp
        p1 = z
        do k = 2, n
          pk = (real(2*k-1,dp)*z*p1 - real(k-1,dp)*p0)/real(k,dp)
          p0 = p1
          p1 = pk
        end do
        pp = real(n,dp)*(z*p1 - p0)/(z*z - 1.0_dp)
        dz = -p1/pp
        z = z + dz
        if (abs(dz) < 1.0e-15_dp) exit
      end do
      x(i) = 0.5_dp*(1.0_dp - z)
      x(n+1-i) = 0.5_dp*(1.0_dp + z)
      wq(i) = 1.0_dp/((1.0_dp - z*z)*pp*pp)
      wq(n+1-i) = wq(i)
    end do

  end subroutine gauss_legendre

  !! X = B A^-1 (solves A^T X^T = B^T row by row)
  function right_div(B, A) result(X)
    real(dp), dimension(nmu,nmu), intent(in) :: B, A
    real(dp), dimension(nmu,nmu) :: X
    real(dp), dimension(nmu,nmu) :: At
    integer, dimension(nmu) :: piv
    integer :: i
    At(:,:) = transpose(A)
    call lu_factor(At, piv)
    do i = 1, nmu
      X(i,:) = B(i,:)
      call lu_solve(At, piv, X(i,:))
    end do
  end function right_div

  !! x = A^-1 b
  function left_solve(A, b) result(x)
    real(dp), dimension(nmu,nmu), intent(in) :: A
    real(dp), dimension(nmu), intent(in) :: b
    real(dp), dimension(nmu) :: x
    real(dp), dimension(nmu,nmu) :: LU
    integer, dimension(nmu) :: piv
    LU(:,:) = A(:,:)
    call lu_factor(LU, piv)
    x(:) = b(:)
    call lu_solve(LU, piv, x)
  end function left_solve

  !! In-place LU factorisation with partial pivoting
  subroutine lu_factor(A, piv)
    real(dp), dimension(nmu,nmu), intent(inout) :: A
    integer, dimension(nmu), intent(out) :: piv
    integer :: i, k, p
    real(dp), dimension(nmu) :: row
    do k = 1, nmu
      p = k - 1 + maxloc(abs(A(k:nmu,k)), dim=1)
      piv(k) = p
      if (p /= k) then
        row(:) = A(k,:)
        A(k,:) = A(p,:)
        A(p,:) = row(:)
      end if
      do i = k+1, nmu
        A(i,k) = A(i,k)/A(k,k)
        A(i,k+1:nmu) = A(i,k+1:nmu) - A(i,k)*A(k,k+1:nmu)
      end do
    end do
  end subroutine lu_factor

  !! Solve LU x = b in place, with the row swaps of lu_factor
  subroutine lu_solve(A, piv, b)
    real(dp), dimension(nmu,nmu), intent(in) :: A
    integer, dimension(nmu), intent(in) :: piv
    real(dp), dimension(nmu), intent(inout) :: b
    integer :: i, k
    real(dp) :: tmp
    do k = 1, nmu
      if (piv(k) /= k) then
        tmp = b(k)
        b(k) = b(piv(k))
        b(piv(k)) = tmp
      end if
    end do
    do i = 2, nmu
      b(i) = b(i) - sum(A(i,1:i-1)*b(1:i-1))
    end do
    do i = nmu, 1, -1
      b(i) = (b(i) - sum(A(i,i+1:nmu)*b(i+1:nmu)))/A(i,i)
    end do
  end subroutine lu_solve

end module sw_AD_mod
