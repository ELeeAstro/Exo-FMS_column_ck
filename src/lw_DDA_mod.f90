!!!
! Elspeth KH Lee - Nov 2026 : Initial version
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!                - Nov 2026 : Optimisations - quadrature tables once per call, even/odd Legendre sums,
!                             sum/difference layer mode outputs, angle independent source function
!                             terms once per layer, thin-layer kernels from the transmissions,
!                             reciprocal angles, contiguous per g-ordinance inputs
!
! lw: Four-stream discrete ordinates adding method (4DDA; Zhang et al. 2013, J. Atmos. Sci.;
!     infrared adding following the Zhang group's four-stream infrared work), thermal emission only.
!     Layer:   analytic 4-stream (2 double-Gauss streams per hemisphere) discrete ordinates
!              solution of each homogeneous layer. With S = I+ + I-, D = I+ - I-, the system
!              reduces to d2S/dtau2 = (alpha+beta)(alpha-beta) S, a 2x2 eigenproblem whose
!              eigenvalues k1^2, k2^2 are the roots of a quadratic (Zhang et al.'s biquadratic).
!              The Planck function is linear in tau inside the layer: I+/- = B(tau) +/- d B1,
!              with (E - w (Ws - Wd)) d = mu, so the thermal particular solution is exact.
!              Only decaying exponentials are used (no overflow for any dtau).
!              The layer is symmetric (same R and T from both sides), so the 4x4 boundary value
!              problem splits into sum and difference parts with E = diag(e^{-k dtau}):
!                R + T = (Gp + Gm E)(Gm + Gp E)^-1,   R - T = (Gp - Gm E)(Gm - Gp E)^-1,
!              and the thermal sources follow from the same two matrices - only 2x2 inverses.
!     Adding:  each layer gives 2x2 reflection R and transmission T (equal from both sides,
!              homogeneous layer) and 2-vector thermal sources (Eu emitted up at the top, Ed
!              down at the bottom). One-sided (Riccati) adding: composites from the bottom
!              (Rbot, Sup, starting from the gas giant boundary I_up = I_down + B(Tint)), then a
!              single downward sweep I_dn(k+1) = (I - R_k Rbot_{k+1})^-1 [T_k I_dn(k) + Ed_k + R_k Sup_{k+1}],
!              I_up(k) = Sup_k + Rbot_k I_dn(k). Identical to two-sided adding, with ~1/3 less work.
!
! Notes:
!   - Henyey-Greenstein phase function with delta-M scaling, f = g**4 (chi_(2 nmu), the first
!     moment the 4-stream expansion cannot hold; the retained moments chi_0..chi_3 are exact).
!   - Double-Gauss quadrature (Gauss-Legendre on [0,1], mu = 1/2 -/+ 1/(2 sqrt(3)), w = 1/2),
!     as Zhang et al.; the streams are the scattering streams, so flux = 2 pi sum(w mu I).
!     With 2 nodes per hemisphere sum_j w_j (P(mu_i,mu_j) + P(mu_i,-mu_j))/2 = 1 exactly,
!     so energy is conserved without renormalising the phase matrix and an isothermal
!     layer gives I = B exactly.
!   - w is capped at w_max (the k2 -> 0 conservative limit is singular in the eigen-solution).
!   - use_sft = .True. (option, off by default): source function technique. The analytic in-layer
!     solution gives S(+/-mu, t) = c0 + c1 t + sum_m [alp_m e^{-k_m t} + bet_m e^{-k_m (D - t)}] at any
!     angle, integrated analytically along the Hogan (2024) Gauss-Jacobi angles (as lw_Fu97 step 2)
!     instead of the 2-node double-Gauss flux integral. Mode coefficients from the level intensities:
!     A + C = M1 [a + b - (B_top + B_bot)], A - C = M2 [a - b - (B_top - B_bot) + 2 d B1].
!     Against 32-stream DISORT (54 layers, Guillot profile, ssa 0-0.99, g 0-0.9): OLR error max
!     0.39% -> 0.15%, worst heating error 12.0 -> 5.4 K/day (ssa = 0: 12.0 -> 2.4 K/day), but heating
!     errors at ssa 0.5-0.9 increase (2.8 -> 5.4 K/day at ssa = 0.75), the 2-node integral cancels
!     part of the four-stream scattering error there. Cost ~1.8x (~1.7x Toon89); more accurate than
!     Toon89 in OLR and heating rate for almost all cases. The 3-node Hogan rule is not worth it
!     (~7% faster, ssa = 0 heating error 2.4 -> 6.1 K/day).
!
! Pros: four-stream accuracy for multiple scattering at two-stream-like cost, no iterations,
!       only 2x2 solves
! Cons: four streams only - for very strong forward scattering lw_AD with more streams is better
!!!

module lw_DDA_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4
  use BB_integrate_mod, only : BB_integrate
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Number of streams per hemisphere (4-stream method) and retained phase function moments
  integer, parameter :: nmu = 2
  integer, parameter :: nmom = 2*nmu

  !! Double-Gauss nodes and weights on [0,1] (sum(w) = 1, sum(w*mu) = 1/2)
  real(dp), dimension(nmu), parameter :: uarr = (/0.2113248654051871_dp, 0.7886751345948129_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.5_dp, 0.5_dp/)
  real(dp), dimension(nmu), parameter :: ruarr = 1.0_dp/uarr

  !! Single scattering albedo cap and minimum layer optical depth
  real(dp), parameter :: w_max = 1.0_dp - 1.0e-6_dp
  real(dp), parameter :: dtau_min = 1.0e-10_dp

  !! Source function technique: flux integration along the Hogan (2024) Gauss-Jacobi-5 angles
  !! (weights for the measure 2 mu dmu, F = pi sum(wh I)) - off by default (canonical 4DDA fluxes)
  logical, parameter :: use_sft = .False.
  integer, parameter :: nh = 4
  real(dp), dimension(nh), parameter :: uh = (/0.0454586727_dp, 0.2322334416_dp, 0.5740198775_dp, 0.9030775973_dp/)
  real(dp), dimension(nh), parameter :: wh = (/0.0092068785_dp, 0.1285704278_dp, 0.4323381850_dp, 0.4298845087_dp/)
  real(dp), dimension(nh), parameter :: ruh = 1.0_dp/uh

  private :: lw_DDA_method, layer_RT, inv2, em1x
  public :: lw_DDA

contains

  subroutine lw_DDA(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
    real(dp), dimension(nlev) :: be_b, tau_g
    real(dp), dimension(nlay) :: ssa_g, gg_g

    !! Quadrature constants (same for every band and g-ordinate)
    integer :: l, i, j, h
    real(dp), dimension(nmom, nmu) :: Plg
    real(dp), dimension(nmom) :: c2l, Ph
    real(dp), dimension(nmom,nmu,nmu) :: PP
    real(dp), dimension(nmom,nmu,nh) :: HPs

    !! Legendre polynomials P_l(mu_i), l = 0..3 (stored at l+1)
    Plg(1,:) = 1.0_dp
    Plg(2,:) = uarr(:)
    Plg(3,:) = 0.5_dp*(3.0_dp*uarr(:)**2 - 1.0_dp)
    Plg(4,:) = 0.5_dp*(5.0_dp*uarr(:)**3 - 3.0_dp*uarr(:))

    !! (2l+1) factors of the phase function expansion
    do l = 0, nmom-1
      c2l(l+1) = real(2*l+1,dp)
    end do

    !! Quadrature-weighted Legendre products, PP(l,i,j) = w_j (2l+1) P_l(mu_i) P_l(mu_j)/2
    do j = 1, nmu
      do i = 1, nmu
        PP(:,i,j) = 0.5_dp*w(j)*c2l(:)*Plg(:,i)*Plg(:,j)
      end do
    end do

    !! Source function technique: the same products with P_l at the Hogan angles
    do h = 1, nh
      Ph(1) = 1.0_dp
      Ph(2) = uh(h)
      Ph(3) = 0.5_dp*(3.0_dp*uh(h)**2 - 1.0_dp)
      Ph(4) = 0.5_dp*(5.0_dp*uh(h)**3 - 3.0_dp*uh(h))
      do j = 1, nmu
        HPs(:,j,h) = 0.5_dp*w(j)*c2l(:)*Ph(:)*Plg(:,j)
      end do
    end do

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
      be_b(:) = be(b,:)
      do g = 1, ng
        tau_g(:) = tau_e(g,b,:)
        ssa_g(:) = ssa(g,b,:)
        gg_g(:) = gg(g,b,:)
        call lw_DDA_method(nlay, nlev, be_b, be_int(b), tau_g, ssa_g, gg_g, PP, HPs, lw_up_g, lw_down_g)
        lw_up(:) = lw_up(:) + lw_up_g(:) * gw(g)
        lw_down(:) = lw_down(:) + lw_down_g(:) * gw(g)
      end do
    end do

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_DDA

  subroutine lw_DDA_method(nlay, nlev, be, be_int, tau_in, w_in, g_in, PP, HPs, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int
    real(dp), dimension(nmom,nmu,nmu), intent(in) :: PP
    real(dp), dimension(nmom,nmu,nh), intent(in) :: HPs

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables and arrays
    integer :: k, l, i, j
    real(dp), dimension(nmom) :: chi
    real(dp), dimension(nlay) :: dtau, w0
    real(dp), dimension(nmu,nmu) :: Ws, Wd, Imat, Mi
    real(dp) :: pe_, po_
    real(dp), dimension(nmu,nmu,nlay) :: R, T, Md
    real(dp), dimension(nmu,nlay) :: Eu, Ed
    real(dp), dimension(nmu,nmu,nlev) :: Rbot
    real(dp), dimension(nmu,nlev) :: Sup
    real(dp), dimension(nmu) :: Iup, Idn

    !! Source function technique work arrays
    integer :: h, m
    real(dp), dimension(nmom, nlay) :: chiL
    real(dp), dimension(nmu, nlay) :: kL, exL, dL
    real(dp), dimension(nmu, nmu, nlay) :: SGL, DGL, M1L, M2L
    real(dp), dimension(nmu, nlev) :: IdnL, IupL
    real(dp), dimension(nmu) :: uu, ss, Ac, Cc, he, ho, alp_p, bet_p, alp_m, bet_m, Kdec, Kgro
    real(dp) :: x, ga, gb, tu, a
    real(dp), dimension(nlev) :: Ihd, Ihu
    real(dp), dimension(nh, nlay) :: Th, up_src, dn_src
    real(dp) :: B1, cb, cd, c0p, c0m, c1, Ti, lin1
    real(dp) :: f

    Imat(:,:) = 0.0_dp
    do i = 1, nmu
      Imat(i,i) = 1.0_dp
    end do

    !! Layer optical depths and single scattering albedo
    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)
    w0(:) = min(max(w_in(:), 0.0_dp), w_max)

    !! Reflection, transmission and thermal sources of each layer
    do k = 1, nlay

      !! Henyey-Greenstein moments, delta-M scaled with the first truncated moment chi_4 = g^4
      do l = 0, nmom-1
        chi(l+1) = g_in(k)**l
      end do
      f = g_in(k)**nmom
      if (w0(k) > 0.0_dp .and. f > 1.0e-12_dp .and. f < 1.0_dp - 1.0e-12_dp) then
        chi(:) = (chi(:) - f)/(1.0_dp - f)
        dtau(k) = (1.0_dp - w0(k)*f)*dtau(k)
        w0(k) = w0(k)*(1.0_dp - f)/(1.0_dp - w0(k)*f)
      end if

      !! Quadrature-weighted phase matrices, normalised so (1/2) int P dmu' = 1:
      !! Ws(i,j) = w_j P(mu_i, mu_j)/2 (same hemisphere), Wd(i,j) = w_j P(mu_i,-mu_j)/2 (opposite)
      !! P_l(-mu) = (-1)^l P_l(mu): Ws = even + odd l terms, Wd = even - odd l terms
      do j = 1, nmu
        do i = 1, nmu
          pe_ = sum(chi(1:nmom:2)*PP(1:nmom:2,i,j))
          po_ = sum(chi(2:nmom:2)*PP(2:nmom:2,i,j))
          Ws(i,j) = pe_ + po_
          Wd(i,j) = pe_ - po_
        end do
      end do

      if (use_sft) chiL(:,k) = chi(:)
      call layer_RT(dtau(k), w0(k), be(k), be(k+1), Ws, Wd, Imat, R(:,:,k), T(:,:,k), Eu(:,k), Ed(:,k), &
        & kL(:,k), exL(:,k), dL(:,k), SGL(:,:,k), DGL(:,:,k), M1L(:,:,k), M2L(:,:,k))

    end do

    !! Upward pass - composites from the bottom: Rbot = reflection of the slab below a level (for
    !! downward light), Sup = upward intensity emitted by the slab below, and Md_k = (I - R_k Rbot_{k+1})^-1
    !! for the downward pass. Lower boundary - gas giant: I_up = I_down + B(Tint) in every stream
    Rbot(:,:,nlev) = Imat(:,:)
    Sup(:,nlev) = be_int
    do k = nlay, 1, -1
      Mi(:,:) = inv2(Imat(:,:) - matmul(Rbot(:,:,k+1), R(:,:,k)))
      Rbot(:,:,k) = R(:,:,k) + matmul(T(:,:,k), matmul(Mi(:,:), matmul(Rbot(:,:,k+1), T(:,:,k))))
      Sup(:,k) = Eu(:,k) + matmul(T(:,:,k), matmul(Mi(:,:), Sup(:,k+1) + matmul(Rbot(:,:,k+1), Ed(:,k))))
      Md(:,:,k) = inv2(Imat(:,:) - matmul(R(:,:,k), Rbot(:,:,k+1)))
    end do

    !! Downward pass - no incident intensity at the top
    Idn(:) = 0.0_dp
    do k = 1, nlev
      if (k > 1) then
        Idn(:) = matmul(Md(:,:,k-1), matmul(T(:,:,k-1), Idn(:)) + Ed(:,k-1) + matmul(R(:,:,k-1), Sup(:,k)))
      end if
      Iup(:) = Sup(:,k) + matmul(Rbot(:,:,k), Idn(:))
      if (use_sft) then
        IdnL(:,k) = Idn(:)
        IupL(:,k) = Iup(:)
      else
        flx_down(k) = 2.0_dp*pi*sum(w(:)*uarr(:)*Idn(:))
        flx_up(k) = 2.0_dp*pi*sum(w(:)*uarr(:)*Iup(:))
      end if
    end do

    if (.not. use_sft) return

    !! Mode coefficients for the actual boundary intensities (I-(0) = IdnL(k), I+(D) = IupL(k+1)):
    !! u = A + C = M1 [a + b - (B_top + B_bot)],  s = A - C = M2 [a - b - (B_top - B_bot) + 2 d B1]
    !! (angle independent - once per layer), then the source terms of each Hogan angle
    do k = 1, nlay
      if (dtau(k) <= dtau_min) then
        Th(:,k) = 1.0_dp
        up_src(:,k) = 0.0_dp
        dn_src(:,k) = 0.0_dp
        cycle
      end if
      B1 = (be(k+1) - be(k))/dtau(k)
      uu(:) = matmul(M1L(:,:,k), IdnL(:,k) + IupL(:,k+1) - (be(k) + be(k+1)))
      ss(:) = matmul(M2L(:,:,k), IdnL(:,k) - IupL(:,k+1) - (be(k) - be(k+1)) + 2.0_dp*dL(:,k)*B1)
      Ac(:) = 0.5_dp*w0(k)*(uu(:) + ss(:))
      Cc(:) = 0.5_dp*w0(k)*(uu(:) - ss(:))
      do h = 1, nh
        !! Phase function projections onto +/-uh, split into even (he) and odd (ho) Legendre terms:
        !! hs = he + ho (same hemisphere), hd = he - ho (opposite)
        do j = 1, nmu
          he(j) = sum(chiL(1:nmom:2,k)*HPs(1:nmom:2,j,h))
          ho(j) = sum(chiL(2:nmom:2,k)*HPs(2:nmom:2,j,h))
        end do
        cb = (1.0_dp - w0(k)) + 2.0_dp*w0(k)*sum(he(:))
        cd = 2.0_dp*w0(k)*sum(ho(:)*dL(:,k))*B1
        c0p = cb*be(k) + cd
        c0m = cb*be(k) - cd
        c1 = cb*B1
        !! Mode amplitudes of the source at +/-uh, with SG = Gp + Gm and DG = Gp - Gm:
        !! hs.Gp + hd.Gm = he.SG + ho.DG,  hs.Gm + hd.Gp = he.SG - ho.DG
        do m = 1, nmu
          ga = sum(he(:)*SGL(:,m,k))
          gb = sum(ho(:)*DGL(:,m,k))
          alp_p(m) = Ac(m)*(ga + gb)
          bet_p(m) = Cc(m)*(ga - gb)
          alp_m(m) = Ac(m)*(ga - gb)
          bet_m(m) = Cc(m)*(ga + gb)
        end do
        tu = dtau(k)*ruh(h)
        Ti = exp(-tu)
        Th(h,k) = Ti
        do m = 1, nmu
          a = kL(m,k)*uh(h)
          Kdec(m) = (1.0_dp - Ti*exL(m,k))/(1.0_dp + a)
          x = a - 1.0_dp
          if (abs(x)*tu > 1.0e-3_dp) then
            Kgro(m) = (Ti - exL(m,k))/x     ! reuses e^{-D/mu} and e^{-k D}, no new exponential
          else
            !! exp(-min(k, 1/uh) dtau) is Ti or e^{-k dtau}; em1x argument |x| tu <= 1e-3 (series)
            Kgro(m) = tu*merge(Ti, exL(m,k), kL(m,k) >= ruh(h)) * em1x(abs(x)*tu)
          end if
        end do
        lin1 = uh(h)*(1.0_dp - Ti) - dtau(k)*Ti
        up_src(h,k) = c0p*(1.0_dp - Ti) + c1*lin1 + sum(alp_p(:)*Kdec(:)) + sum(bet_p(:)*Kgro(:))
        dn_src(h,k) = (c0m + c1*dtau(k))*(1.0_dp - Ti) - c1*lin1 + sum(alp_m(:)*Kgro(:)) + sum(bet_m(:)*Kdec(:))
      end do
    end do

    !! Formal solution along each Hogan angle
    flx_down(:) = 0.0_dp
    flx_up(:) = 0.0_dp
    do h = 1, nh
      Ihd(1) = 0.0_dp
      do k = 1, nlay
        Ihd(k+1) = Ihd(k)*Th(h,k) + dn_src(h,k)
      end do
      Ihu(nlev) = Ihd(nlev) + be_int
      do k = nlay, 1, -1
        Ihu(k) = Ihu(k+1)*Th(h,k) + up_src(h,k)
      end do
      flx_down(:) = flx_down(:) + pi*wh(h)*Ihd(:)
      flx_up(:) = flx_up(:) + pi*wh(h)*Ihu(:)
    end do

  end subroutine lw_DDA_method

  !! Analytic 4-stream discrete ordinates solution of one homogeneous layer:
  !! R, T (2x2, incident from either side) and thermal sources Eu (up, top), Ed (down, bottom)
  subroutine layer_RT(dtau, w0, B_top, B_bot, Ws, Wd, Imat, R, T, Eu, Ed, kko, exo, dout, SGo, DGo, M1o, M2o)
    implicit none

    real(dp), intent(in) :: dtau, w0, B_top, B_bot
    real(dp), dimension(nmu,nmu), intent(in) :: Ws, Wd, Imat
    real(dp), dimension(nmu,nmu), intent(out) :: R, T
    real(dp), dimension(nmu), intent(out) :: Eu, Ed
    real(dp), dimension(nmu), intent(out) :: kko, exo, dout
    real(dp), dimension(nmu,nmu), intent(out) :: SGo, DGo, M1o, M2o

    integer :: i, m
    real(dp), dimension(nmu,nmu) :: al, bt, amb, apb, Q, Gp, Gm, GmE, GpE, P1, P2
    real(dp), dimension(nmu) :: kk, ex, X, Y, v1, v2, d, one, sum_o, dif_o
    real(dp) :: tr, det, disc, lam1, lam, B1

    !! Optically negligible layer - transparent
    if (dtau <= dtau_min) then
      R(:,:) = 0.0_dp
      T(:,:) = Imat(:,:)
      Eu(:) = 0.0_dp
      Ed(:) = 0.0_dp
      kko(:) = 0.0_dp ; exo(:) = 1.0_dp ; dout(:) = 0.0_dp
      SGo(:,:) = 0.0_dp ; DGo(:,:) = 0.0_dp ; M1o(:,:) = 0.0_dp ; M2o(:,:) = 0.0_dp
      return
    end if

    !! mu dI+/dtau = alpha I+ - beta I-,  dI-/dtau = beta I+ - alpha I-  (per unit tau, divided by mu)
    do i = 1, nmu
      al(i,:) = (Imat(i,:) - w0*Ws(i,:))*ruarr(i)
      bt(i,:) = w0*Wd(i,:)*ruarr(i)
    end do
    amb(:,:) = al(:,:) - bt(:,:)
    apb(:,:) = al(:,:) + bt(:,:)
    Q(:,:) = matmul(apb(:,:), amb(:,:))

    !! Eigenvalues k^2 of Q (quadratic in k^2); smaller root via det/lam1 to avoid cancellation
    tr = Q(1,1) + Q(2,2)
    det = Q(1,1)*Q(2,2) - Q(1,2)*Q(2,1)
    disc = sqrt(max(tr**2 - 4.0_dp*det, 0.0_dp))
    lam1 = 0.5_dp*(tr + disc)
    kk(1) = sqrt(lam1)
    kk(2) = sqrt(max(det/lam1, tiny(1.0_dp)))

    !! Eigenvectors X of Q, then the decaying (e^{-k tau}) modes:
    !! I+ = Gp e^{-k tau}, I- = Gm e^{-k tau}, with Gp = (X - Y)/2, Gm = (X + Y)/2, Y = (alpha-beta) X/k
    !! The growing mode is the mirror image (I+ = Gm, I- = Gp), written as e^{-k (dtau - tau)}
    do m = 1, nmu
      lam = kk(m)**2
      v1(:) = (/Q(1,2), lam - Q(1,1)/)
      v2(:) = (/lam - Q(2,2), Q(2,1)/)
      if (sum(v1**2) >= sum(v2**2)) then
        X(:) = v1(:)*(1.0_dp/sqrt(sum(v1**2)))
      else
        X(:) = v2(:)*(1.0_dp/sqrt(sum(v2**2)))
      end if
      Y(:) = matmul(amb(:,:), X(:))*(1.0_dp/kk(m))
      Gp(:,m) = 0.5_dp*(X(:) - Y(:))
      Gm(:,m) = 0.5_dp*(X(:) + Y(:))
      ex(m) = exp(-kk(m)*dtau)
    end do

    !! Thermal particular solution, B(tau) = B_top + B1 tau: I+/- = B(tau) +/- d B1
    B1 = (B_bot - B_top)/dtau
    d(:) = matmul(inv2(Imat(:,:) - w0*(Ws(:,:) - Wd(:,:))), uarr(:))

    !! Boundary value problem for the coefficients (A_m, C_m):
    !! I+ = sum_m [A_m Gp_m e^{-k_m tau} + C_m Gm_m e^{-k_m (dtau - tau)}] + Ip+
    !! I- = sum_m [A_m Gm_m e^{-k_m tau} + C_m Gp_m e^{-k_m (dtau - tau)}] + Ip-
    !! with I-(0) = a and I+(dtau) = b given. With u = A + C, s = A - C:
    !!   a + b - (Ip-(0) + Ip+(D)) = (Gm + Gp E) u,   a - b - (Ip-(0) - Ip+(D)) = (Gm - Gp E) s
    !!   I+(0) + I-(D) = (Gp + Gm E) u + Ip+(0) + Ip-(D),   I+(0) - I-(D) = (Gp - Gm E) s + Ip+(0) - Ip-(D)
    do m = 1, nmu
      GmE(:,m) = Gm(:,m)*ex(m)
      GpE(:,m) = Gp(:,m)*ex(m)
    end do
    M1o(:,:) = inv2(Gm(:,:) + GpE(:,:))
    M2o(:,:) = inv2(Gm(:,:) - GpE(:,:))
    P1(:,:) = matmul(Gp(:,:) + GmE(:,:), M1o(:,:))   ! R + T
    P2(:,:) = matmul(Gp(:,:) - GmE(:,:), M2o(:,:))   ! R - T
    if (use_sft) then
      kko(:) = kk(:) ; exo(:) = ex(:) ; dout(:) = d(:) ; SGo(:,:) = Gp(:,:) + Gm(:,:) ; DGo(:,:) = Gp(:,:) - Gm(:,:)
    end if

    R(:,:) = 0.5_dp*(P1(:,:) + P2(:,:))
    T(:,:) = 0.5_dp*(P1(:,:) - P2(:,:))

    !! Thermal sources (a = b = 0): Ip-(0) + Ip+(D) = Ip+(0) + Ip-(D) = (B_top + B_bot),
    !! Ip-(0) - Ip+(D) = (B_top - B_bot) - 2 d B1,  Ip+(0) - Ip-(D) = (B_top - B_bot) + 2 d B1
    one(:) = 1.0_dp
    sum_o(:) = (B_top + B_bot)*(one(:) - matmul(P1(:,:), one(:)))
    dif_o(:) = (B_top - B_bot)*(one(:) - matmul(P2(:,:), one(:))) + 2.0_dp*B1*(d(:) + matmul(P2(:,:), d(:)))
    Eu(:) = 0.5_dp*(sum_o(:) + dif_o(:))
    Ed(:) = 0.5_dp*(sum_o(:) - dif_o(:))

  end subroutine layer_RT

  !! (1 - exp(-x))/x, stable as x -> 0
  elemental function em1x(x) result(r)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: r
    if (x < 1.0e-3_dp) then
      r = 1.0_dp - x*(0.5_dp - x*(1.0_dp/6.0_dp - x/24.0_dp))
    else
      r = (1.0_dp - exp(-x))/x
    end if
  end function em1x

  !! Inverse of a 2x2 matrix
  pure function inv2(A) result(B)
    implicit none

    real(dp), dimension(2,2), intent(in) :: A
    real(dp), dimension(2,2) :: B
    real(dp) :: rdet

    rdet = 1.0_dp/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
    B(1,1) = A(2,2)*rdet
    B(1,2) = -A(1,2)*rdet
    B(2,1) = -A(2,1)*rdet
    B(2,2) = A(1,1)*rdet

  end function inv2

end module lw_DDA_mod