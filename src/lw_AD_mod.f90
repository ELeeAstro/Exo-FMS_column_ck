!!!
! Elspeth KH Lee - Jun 2024 : Initial version
!                - Nov 2026 : Full implementation, converted from AddingDoublingRT (JAX/C++)
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!
! lw: Adding-doubling method (matrix operator method; Plass, Kattawar & Catchings 1973;
!     Hansen 1971), converted to Fortran from AddingDoublingRT (src_jax/doubling.py,
!     adding.py, solver.py; same algorithm as the C++ src/), thermal emission only.
!     Doubling:  reflection R, transmission T and thermal source vectors of each homogeneous
!                layer, started from an optically thin sub-layer (single scattering) and
!                doubled nn times, with the Planck function linear in tau inside the layer.
!     Adding:    composites from the top (rtop) and from the bottom (rbase) of the column,
!                then the up and down intensities at every interface from the two composites.
!
! Pros: accurate multiple scattering for any w and g, converges with nmu
! Cons: nn ~ log2(dtau) + 16 doublings (N x N matrix solves) per layer - the most expensive
!       of the matrix-free/two-stream class, though far cheaper than lw_disort
!!!

module lw_AD_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4
  use BB_integrate_mod, only : BB_integrate
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Number of Gauss-Legendre streams per hemisphere (2*nmu stream method)
  integer, parameter :: nmu = 4

  !! Number of phase function Legendre moments retained (delta-M truncates chi_(2 nmu))
  integer, parameter :: nmom = 2*nmu

  !! Layer (or composite) reflection, transmission and thermal sources.
  !! ab: incident from the top (a) side, ba: incident from the bottom (b) side.
  type :: ad_layer
    real(dp), dimension(nmu,nmu) :: Rab, Rba, Tab, Tba
    real(dp), dimension(nmu) :: sup, sdn
  end type ad_layer

  private :: lw_adding_doubling, layer_doubling, add_layers, gauss_legendre, &
    & right_div, left_solve, lu_factor, lu_solve, ipow0
  public :: lw_AD

contains

  subroutine lw_AD(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev, nb, ng
    real(dp), dimension(ng), intent(in) :: gw
    real(dp), dimension(nb+1), intent(in) :: wn_e
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(ng,nb,nlev), intent(in) :: tau_e
    real(dp), dimension(ng,nb,nlay), intent(in) :: ssa, gg
    real(dp), intent(in) :: Tint

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: lw_up, lw_down, lw_net
    real(dp), intent(out) :: olr

    !! Work variables
    integer :: k, b, g
    real(dp), dimension(nlev) :: Te
    real(dp), dimension(nb,nlev) :: be
    real(dp), dimension(nb) :: be_int
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    !! Use WENO4 method to (smoothly) interpolate layers to levels, in ln p
    Te(:) = interpolate_weno4(log(pe), log(pl), Tl, .False.)

    !! Edges are linearly extrapolated in log(p)-log(T) space.
    Te(1) = 10.0_dp**(log10(Te(2)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * &
      & log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Te(nlay)) + (log10(pe(nlev)/pe(nlay))/ &
      & log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    !! Find integrated planck function for each band for each level
    do k = 1, nlev
      call BB_integrate(nb, Te(k), wn_e, be(:,k))
    end do
    call BB_integrate(nb, Tint, wn_e, be_int) ! Integrated planck function intensity for internal temperature

    !! Longwave flux calculation - sum over each band and g-ordinance
    lw_up(:) = 0.0_dp
    lw_down(:) = 0.0_dp
    do b = 1, nb
      do g = 1, ng
        call lw_adding_doubling(nlay, nlev, be(b,:), be_int(b), tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), &
          & lw_up_g(:), lw_down_g(:))
        lw_up(:) = lw_up(:) + lw_up_g(:) * gw(g)
        lw_down(:) = lw_down(:) + lw_down_g(:) * gw(g)
      end do
    end do

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_AD

  subroutine lw_adding_doubling(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables and arrays
    integer :: k, l, i, j, ltot
    real(dp), dimension(nmu) :: mu, wt, Iup, Idn
    real(dp), dimension(nmom, nmu) :: Pl
    real(dp), dimension(nmom) :: chi, c2l, sgn
    real(dp), dimension(nmu,nmu) :: Ppp, Ppm, Imat
    real(dp), dimension(nlay) :: dtau, w0
    real(dp) :: f, colsum
    type(ad_layer), dimension(nlay+1) :: lay
    type(ad_layer), dimension(nlay+1) :: rtop, rbase

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

    !! Layer optical depths
    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)
    w0(:) = min(max(w_in(:), 0.0_dp), 1.0_dp)

    !! Doubling - R, T and thermal sources of each layer
    do k = 1, nlay

      Ppp(:,:) = 0.0_dp
      Ppm(:,:) = 0.0_dp

      if (w0(k) > 0.0_dp .and. dtau(k) > 0.0_dp) then

        !! Henyey-Greenstein moments, delta-M scaled with the first truncated moment chi_(2 nmu)
        f = g_in(k)**nmom
        do l = 0, nmom-1
          chi(l+1) = g_in(k)**l
        end do
        if (f > 1.0e-12_dp .and. f < 1.0_dp - 1.0e-12_dp) then
          chi(:) = (chi(:) - f)/(1.0_dp - f)
          dtau(k) = (1.0_dp - w0(k)*f)*dtau(k)
          w0(k) = w0(k)*(1.0_dp - f)/(1.0_dp - w0(k)*f)
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

      end if

      call layer_doubling(dtau(k), w0(k), be(k), be(k+1), Ppp, Ppm, mu, wt, Imat, lay(k))

    end do

    !! Lower boundary "surface layer": I_up = I_down + B(Tint) in every stream
    ltot = nlay + 1
    lay(ltot)%Rab(:,:) = Imat(:,:)
    lay(ltot)%Rba(:,:) = Imat(:,:)
    lay(ltot)%Tab(:,:) = 0.0_dp
    lay(ltot)%Tba(:,:) = 0.0_dp
    lay(ltot)%sup(:) = be_int
    lay(ltot)%sdn(:) = 0.0_dp

    !! Adding - composites from the bottom (rbase(n): bottom n layers) and from the top
    !! (rtop(n): top n layers)
    rbase(1) = lay(ltot)
    do l = 2, ltot
      call add_layers(lay(ltot-l+1), rbase(l-1), Imat, rbase(l))
    end do

    rtop(1) = lay(1)
    do l = 2, nlay
      call add_layers(rtop(l-1), lay(l), Imat, rtop(l))
    end do

    !! Top of atmosphere - no incident intensity from above
    Iup(:) = rbase(ltot)%sup(:)
    flx_up(1) = 2.0_dp*pi*sum(wt(:)*mu(:)*Iup(:))
    flx_down(1) = 0.0_dp

    !! Interfaces below layer l: top composite rtop(l), base composite rbase(ltot-l)
    !! I_up = (I - Rab_b Rba_t)^-1 (Rab_b sdn_t + sup_b),  I_dn = (I - Rba_t Rab_b)^-1 (Rba_t sup_b + sdn_t)
    do l = 1, nlay
      associate (t => rtop(l), b => rbase(ltot-l))
        Iup(:) = left_solve(Imat - matmul(b%Rab, t%Rba), matmul(b%Rab, t%sdn) + b%sup)
        Idn(:) = left_solve(Imat - matmul(t%Rba, b%Rab), matmul(t%Rba, b%sup) + t%sdn)
      end associate
      flx_up(l+1) = 2.0_dp*pi*sum(wt(:)*mu(:)*Iup(:))
      flx_down(l+1) = 2.0_dp*pi*sum(wt(:)*mu(:)*Idn(:))
    end do

  end subroutine lw_adding_doubling

  !! R, T and thermal sources of one homogeneous layer by doubling
  subroutine layer_doubling(tau, omega, B_top, B_bot, Ppp, Ppm, mu, wt, Imat, lyr)
    implicit none

    real(dp), intent(in) :: tau, omega, B_top, B_bot
    real(dp), dimension(nmu,nmu), intent(in) :: Ppp, Ppm, Imat
    real(dp), dimension(nmu), intent(in) :: mu, wt
    type(ad_layer), intent(out) :: lyr

    integer :: i, n, nn
    real(dp) :: B_bar, B_d, tau0, g, con
    real(dp), dimension(nmu) :: trans, slope, y, z, zpgy
    real(dp), dimension(nmu,nmu) :: R, T, TG, TGR, Gpp, Gpm

    !! Planck function linear in tau: mean and slope over the layer
    B_bar = 0.5_dp*(B_bot + B_top)
    if (tau > 0.0_dp) then
      B_d = (B_bot - B_top)/tau
    else
      B_d = 0.0_dp
    end if

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

    !! Pure absorption - exact diagonal transmission and linear-in-tau emission
    if (omega <= 0.0_dp) then
      trans(:) = exp(max(-tau/mu(:), -200.0_dp))
      slope(:) = mu(:)*(1.0_dp - trans(:)) - 0.5_dp*tau*(1.0_dp + trans(:))
      lyr%Tab(:,:) = 0.0_dp
      do i = 1, nmu
        lyr%Tab(i,i) = trans(i)
      end do
      lyr%Tba(:,:) = lyr%Tab(:,:)
      lyr%sup(:) = B_bar*(1.0_dp - trans(:)) + B_d*slope(:)
      lyr%sdn(:) = B_bar*(1.0_dp - trans(:)) - B_d*slope(:)
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

    !! Initial thin sub-layer: R, T and the thermal source basis (s = y B_bar +/- z B_d)
    R(:,:) = tau0*Gpm(:,:)
    T(:,:) = Imat(:,:) - tau0*Gpp(:,:)
    y(:) = (1.0_dp - omega)*tau0/mu(:)
    z(:) = 0.0_dp
    g = 0.5_dp*tau0

    !! Doubling: combine two identical sub-layers of depth 2g into one of depth 4g
    do n = 1, nn
      TG(:,:) = right_div(T, Imat - matmul(R, R))  ! T (I - R^2)^-1
      TGR(:,:) = matmul(TG, R)

      zpgy(:) = z(:) + g*y(:)
      z(:) = matmul(TG - TGR, zpgy) + z(:) - g*y(:)
      y(:) = matmul(TG + TGR, y) + y(:)

      R(:,:) = R(:,:) + matmul(TGR, T)
      T(:,:) = matmul(TG, T)

      g = 2.0_dp*g
    end do

    !! Homogeneous layer: symmetric about its mid-plane
    lyr%Rab(:,:) = R(:,:)
    lyr%Rba(:,:) = R(:,:)
    lyr%Tab(:,:) = T(:,:)
    lyr%Tba(:,:) = T(:,:)
    lyr%sup(:) = y(:)*B_bar + z(:)*B_d
    lyr%sdn(:) = y(:)*B_bar - z(:)*B_d

  end subroutine layer_doubling

  !! Combine a top layer and a bottom layer into one composite
  subroutine add_layers(top, bot, Imat, comb)
    implicit none

    type(ad_layer), intent(in) :: top, bot
    real(dp), dimension(nmu,nmu), intent(in) :: Imat
    type(ad_layer), intent(out) :: comb

    real(dp), dimension(nmu,nmu) :: TD1, TD2

    !! Multiple reflections between the two layers: (I - Rab_bot Rba_top)^-1, (I - Rba_top Rab_bot)^-1
    TD1(:,:) = right_div(top%Tba, Imat - matmul(bot%Rab, top%Rba))
    TD2(:,:) = right_div(bot%Tab, Imat - matmul(top%Rba, bot%Rab))

    comb%Rab(:,:) = top%Rab(:,:) + matmul(matmul(TD1, bot%Rab), top%Tab)
    comb%Rba(:,:) = bot%Rba(:,:) + matmul(matmul(TD2, top%Rba), bot%Tba)
    comb%Tab(:,:) = matmul(TD2, top%Tab)
    comb%Tba(:,:) = matmul(TD1, bot%Tba)

    comb%sup(:) = top%sup(:) + matmul(TD1, bot%sup + matmul(bot%Rab, top%sdn))
    comb%sdn(:) = bot%sdn(:) + matmul(TD2, top%sdn + matmul(top%Rba, bot%sup))

  end subroutine add_layers

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

end module lw_AD_mod
