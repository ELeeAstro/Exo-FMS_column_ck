!!!
! Elspeth KH Lee - Nov 2026 : Initial version
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!                - Nov 2026 : Optimisations - Legendre tables once per call, layer-last array layout,
!                             direct beam as a running product of the layer transmissions,
!                             reciprocals, contiguous per g-ordinance inputs
!
! sw: delta-four-stream discrete ordinates method (Liou, Fu & Ackerman 1988, J. Atmos. Sci.)
!     Four streams at the double-Gauss points (mu1, mu2) = (0.2113, 0.7887), a1 = a2 = 1/2,
!     four-term phase function expansion (l = 0..3) with delta-M scaling (f = chi_4, g**4 for HG).
!     Each layer is solved analytically:
!       - the 4x4 eigenproblem reduces (symmetric quadrature) to the 2x2 problem C v = k**2 v,
!         C = (A - B)(A + B), so k1**2 and k2**2 are the roots of a quadratic (Liou et al. Eq. 14)
!       - the direct beam particular solution reduces to the 2x2 system
!         (C - I/mu0**2) Z_s = (A - B) M^-1 Q_s - M^-1 Q_d / mu0
!     The layer reflection and transmission (2x2, in stream intensities) and the direct beam
!     source terms follow from the sum/difference split of the boundary conditions (2x2 inverses only).
!     Layers are combined with one-sided (Riccati) adding (gas giant lower boundary, no reflection):
!     an upward pass for the reflection of everything below each level (Rb) and the beam source (RkN),
!     then a single downward sweep D_{l+1} = (I - R_l Rb_{l+1})^-1 [T_l D_l + Tdir_l T_l + R_l RkN_{l+1} T_{l+1}],
!     U_k = Rb_k D_k + RkN_k T_k. Identical to two-sided adding with ~1/3 less work; no banded matrix.
!
! Notes:
!   - Liou et al. (1988) solve the whole column with a 4N banded matrix (as sw_SH4); the adding
!     here gives the same answer for the same layer solutions.
!   - Energy is conserved exactly: the double-Gauss rule integrates P_l(mu) P_l(mu') (l <= 3) over
!     each hemisphere exactly.
!   - w is capped at w_max = 1 - 1e-12 (k -> 0 conservative limit), the eigenvector difference vector is
!     found from (A - B)^-1 so that the near-conservative limit stays well conditioned. (A 1 - 1e-6 cap
!     absorbs ~0.2% of the beam at w = 1 in a tau ~ 6000 column; 1 - 1e-12 agrees with sw_AD to 1e-5.)
!   - mu0 is nudged off the k*mu0 = 1 resonance as in sw_SDA and sw_Toon89.
!
! Pros: Four-stream accuracy with closed-form eigenvalues and only 2x2 matrix inversions
! Cons: Four-stream approximation
!!!

module sw_Liou88_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: fourpi = 4.0_dp * pi
  real(dp), parameter :: r_fourpi = 1.0_dp/fourpi

  !! Double-Gauss streams and weights on [0,1] (sum(a) = 1, flux = 2 pi sum(a mu I))
  integer, parameter :: nmu = 2
  real(dp), dimension(nmu), parameter :: uarr = (/0.2113248654051871_dp, 0.7886751345948129_dp/)
  real(dp), dimension(nmu), parameter :: wa = (/0.5_dp, 0.5_dp/)
  real(dp), dimension(nmu), parameter :: ruarr = 1.0_dp/uarr

  !! Number of phase function moments retained (chi_0..chi_3), delta-M uses chi_4
  integer, parameter :: nstr = 4

  !! Single scattering albedo cap
  real(dp), parameter :: w_max = 1.0_dp - 1.0e-12_dp

  !! How close to the k*mu0 = 1 resonance before mu0 is nudged
  real(dp), parameter :: res_eps = 1.0e-6_dp

  private :: sw_Liou88_method, eigvec2, matinv2, Pleg
  public :: sw_Liou88

contains

  subroutine sw_Liou88(nlay, nlev, nb, ng, gw, tau_e, mu_z, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev, nb, ng
    real(dp), dimension(ng), intent(in) :: gw
    real(dp), dimension(ng,nb,nlev), intent(in) :: tau_e
    real(dp), dimension(ng,nb,nlay), intent(in) :: ssa, gg
    real(dp), intent(in) :: mu_z
    real(dp), dimension(nb), intent(in) :: Finc

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down, sw_net
    real(dp), intent(out) :: asr

    !! Work variables
    integer :: b, g, i, j, l
    real(dp), dimension(nlev) :: sw_up_g, sw_down_g
    real(dp), dimension(nlev) :: tau_g
    real(dp), dimension(nlay) :: ssa_g, gg_g
    real(dp), dimension(0:nstr-1,nmu) :: PlU
    real(dp), dimension(0:nstr-1) :: Pl0, sgnl
    real(dp), dimension(0:nstr-1,nmu,nmu) :: PS

    sw_up(:) = 0.0_dp
    sw_down(:) = 0.0_dp

    !! Shortwave flux calculation - sum over each band and g-ordinance
    if (mu_z > 0.0_dp) then

      !! Legendre polynomials P_l(mu_i) at the streams and P_l(mu0), (-1)^l, and
      !! PS(l,i,j) = a_j P_l(mu_i) P_l(mu_j) - the same for every band and g-ordinance
      do i = 1, nmu
        do l = 0, nstr-1
          PlU(l,i) = Pleg(l, uarr(i))
        end do
      end do
      do l = 0, nstr-1
        sgnl(l) = real(1 - 2*mod(l,2), dp)
        Pl0(l) = Pleg(l, mu_z)
      end do
      do j = 1, nmu
        do i = 1, nmu
          PS(:,i,j) = wa(j)*PlU(:,i)*PlU(:,j)
        end do
      end do

      do b = 1, nb
        if (Finc(b) <= 0.0_dp) cycle
        do g = 1, ng
          ! Contiguous copies (avoids array temporaries for the strided slices)
          tau_g(:) = tau_e(g,b,:)
          ssa_g(:) = ssa(g,b,:)
          gg_g(:) = gg(g,b,:)
          call sw_Liou88_method(nlay, nlev, Finc(b), mu_z, tau_g, ssa_g, gg_g, PlU, Pl0, sgnl, PS, &
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

  end subroutine sw_Liou88

  subroutine sw_Liou88_method(nlay, nlev, F0_in, mu_in, tau_in, w_in, g_in, PlU, Pl0, sgnl, PS, flx_down, flx_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0_in, mu_in
    real(dp), dimension(nlev), intent(in) :: tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Legendre tables (from the wrapper): P_l(mu_i), P_l(mu0), (-1)^l, a_j P_l(mu_i) P_l(mu_j)
    real(dp), dimension(0:nstr-1,nmu), intent(in) :: PlU
    real(dp), dimension(0:nstr-1), intent(in) :: Pl0, sgnl
    real(dp), dimension(0:nstr-1,nmu,nmu), intent(in) :: PS

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    !! Work variables
    integer :: k, l, i, j, m, lp1
    real(dp) :: mu_z, mu_l, f0, trC, detC, disc, rfc, w4
    real(dp), dimension(nlay) :: dtau, w0, fc, dtr, dtr_z
    real(dp), dimension(nlev) :: T
    real(dp), dimension(0:nstr) :: chi
    real(dp), dimension(0:nstr-1) :: om
    real(dp), dimension(nmu) :: kk, ee, v, d, Qp, Qm, Qs, Qd, Zs, Zd, Zp, Zm, rhs
    real(dp), dimension(nmu) :: sum_o, dif_o
    real(dp), dimension(nmu,nmu) :: ApB, AmB, AmB_i, C, X, Y, XE, YE, E, Wk

    real(dp), dimension(nmu,nlay) :: Rdir, Tdir
    real(dp), dimension(nmu,nmu,nlay) :: Rdiff, Tdiff

    real(dp), dimension(nmu,nlev) :: RkN, U_I, D_I
    real(dp), dimension(nmu,nmu,nlev) :: RbkN
    real(dp), dimension(nmu,nmu,nlay) :: Minv
    real(dp), dimension(nmu,nmu) :: TT, CC
    real(dp), dimension(nmu) :: DD

    real(dp), dimension(0:nstr-1) :: Pl0l
    real(dp) :: se, so

    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w_in(:) <= 1.0e-12_dp)) then

      flx_down(:) = F0_in*mu_in * exp(-(tau_in(:) - tau_in(1))/mu_in)

      flx_up(:) = 0.0_dp

      return

    end if

    !! Identity matrix
    E(:,:) = 0.0_dp
    E(1,1) = 1.0_dp ; E(2,2) = 1.0_dp

    !! Calculate dtau in each layer
    dtau(:) = tau_in(2:) - tau_in(1:nlay)

    !! Single zenith angle for the whole column (value at the lower boundary)
    mu_z = mu_in

    !! Delta-M scaling (f = chi_4 = g**4 for HG, as Liou et al. 1988 and sw_SDA)
    fc(:) = g_in(:)**nstr
    w0(:) = min(w_in(:)*((1.0_dp - fc(:))/(1.0_dp - fc(:)*w_in(:))), w_max)
    dtau(:) = (1.0_dp - w_in(:)*fc(:))*dtau(:)

    !! Find the reflection and transmission coefficients (direct and diffuse) for each layer
    do k = 1, nlay

      !! Scaled phase function expansion coefficients, omega_l = (2l+1)(chi_l - f)/(1 - f)
      chi(0) = 1.0_dp
      do l = 1, nstr
        chi(l) = chi(l-1)*g_in(k)
      end do
      rfc = 1.0_dp/(1.0_dp - fc(k))
      do l = 0, nstr-1
        om(l) = real(2*l+1,dp) * (chi(l) - fc(k))*rfc
      end do

      !! Scattering matrices Ppp(i,j) = w/2 a_j P(mu_i,mu_j), Ppm(i,j) = w/2 a_j P(mu_i,-mu_j):
      !! Ppp + Ppm = w sum_even omega_l PS_l,  Ppp - Ppm = w sum_odd omega_l PS_l
      !! A +/- B, with A = M^-1 (I - Ppp), B = -M^-1 Ppm
      do j = 1, nmu
        do i = 1, nmu
          se = om(0)*PS(0,i,j) + om(2)*PS(2,i,j)
          so = om(1)*PS(1,i,j) + om(3)*PS(3,i,j)
          ApB(i,j) = (E(i,j) - w0(k)*se)*ruarr(i)
          AmB(i,j) = (E(i,j) - w0(k)*so)*ruarr(i)
        end do
      end do
      AmB_i(:,:) = matinv2(AmB(:,:))
      C(:,:) = matmul(AmB(:,:), ApB(:,:))

      !! Eigenvalues k**2 of the 2x2 C matrix - closed form quadratic
      trC = C(1,1) + C(2,2)
      detC = C(1,1)*C(2,2) - C(1,2)*C(2,1)
      disc = sqrt(max(trC**2 - 4.0_dp*detC, 0.0_dp))
      kk(1) = sqrt(max(0.5_dp*(trC + disc), 0.0_dp))
      kk(2) = sqrt(max(detC/(0.5_dp*(trC + disc)), 0.0_dp)) ! smaller root without cancellation

      !! Homogeneous modes I+ = X e^{-k t}, I- = Y e^{-k t} (the e^{-k (D-t)} mode swaps X and Y)
      !! S = I+ + I- = v, D = I+ - I- = d = -k (A - B)^-1 v
      do m = 1, nmu
        call eigvec2(C(:,:), kk(m)**2, m, v(:))
        d(:) = -kk(m)*matmul(AmB_i(:,:), v(:))
        X(:,m) = 0.5_dp*(v(:) + d(:))
        Y(:,m) = 0.5_dp*(v(:) - d(:))
        ee(m) = exp(-kk(m)*dtau(k))
      end do

      !! Zenith angle for this layer - nudge mu0 off the k*mu0 = 1 resonance (0/0 in the particular solution)
      mu_l = mu_z
      if (abs(1.0_dp - kk(1)*mu_l) < res_eps) mu_l = mu_l + 2.0_dp*res_eps/kk(1)
      if (abs(1.0_dp - kk(2)*mu_l) < res_eps) mu_l = mu_l + 2.0_dp*res_eps/kk(2)
      f0 = 1.0_dp/mu_l

      !! Layer direct beam transmission (nudged mu0), and with the column mu0 for the beam itself
      dtr(k) = exp(-dtau(k)*f0)
      if (mu_l /= mu_z) then
        dtr_z(k) = exp(-dtau(k)/mu_z)
      else
        dtr_z(k) = dtr(k)
      end if

      !! Direct beam source - normalise the beam flux (normal to the beam) to 1 here
      !! Q+_i = w/(4 pi) P(mu_i,-mu0) = w/(4 pi) sum (-1)^l omega_l P_l(mu_i) P_l(mu0), Q-_i: no (-1)^l
      if (mu_l /= mu_z) then
        do l = 0, nstr-1
          Pl0l(l) = Pleg(l, mu_l)
        end do
      else
        Pl0l(:) = Pl0(:)
      end if
      w4 = w0(k)*r_fourpi
      do i = 1, nmu
        Qp(i) = w4 * sum(sgnl(:)*om(:)*PlU(:,i)*Pl0l(:))
        Qm(i) = w4 * sum(om(:)*PlU(:,i)*Pl0l(:))
      end do
      Qs(:) = (Qp(:) + Qm(:))*ruarr(:)
      Qd(:) = (Qp(:) - Qm(:))*ruarr(:)

      !! Particular solution I+/- = Z+/- e^{-t/mu0}
      Wk(:,:) = C(:,:) - E(:,:)*f0**2
      rhs(:) = matmul(AmB(:,:), Qs(:)) - f0*Qd(:)
      Zs(:) = matmul(matinv2(Wk(:,:)), rhs(:))
      Zd(:) = -(matmul(ApB(:,:), Zs(:)) - Qs(:))*mu_l
      Zp(:) = 0.5_dp*(Zs(:) + Zd(:))
      Zm(:) = 0.5_dp*(Zs(:) - Zd(:))

      !! Boundary conditions I-(0) = I_top, I+(D) = I_bot, split into sum and difference:
      !!  I_top + I_bot = (Y + X E)(alpha + beta),  I_top - I_bot = (Y - X E)(alpha - beta)
      !!  I+(0) + I-(D) = (X + Y E)(alpha + beta),  I+(0) - I-(D) = (X - Y E)(alpha - beta)
      do m = 1, nmu
        XE(:,m) = X(:,m)*ee(m)
        YE(:,m) = Y(:,m)*ee(m)
      end do
      TT(:,:) = matmul(X(:,:) + YE(:,:), matinv2(Y(:,:) + XE(:,:)))  ! R + T
      CC(:,:) = matmul(X(:,:) - YE(:,:), matinv2(Y(:,:) - XE(:,:)))  ! R - T

      !! Diffuse reflection and transmission (intensity, identical from both sides)
      Rdiff(:,:,k) = 0.5_dp*(TT(:,:) + CC(:,:))
      Tdiff(:,:,k) = 0.5_dp*(TT(:,:) - CC(:,:))

      !! Direct beam: no incident diffuse intensity, I-(0) = 0 and I+(D) = 0, so
      !! alpha + beta = -(Y + X E)^-1 (Z- + Z+ T0) and the outputs reuse R + T and R - T
      sum_o(:) = Zp(:) + Zm(:)*dtr(k) - matmul(TT(:,:), Zm(:) + Zp(:)*dtr(k))
      dif_o(:) = Zp(:) - Zm(:)*dtr(k) - matmul(CC(:,:), Zm(:) - Zp(:)*dtr(k))

      !! Direct beam reflection and transmission per unit beam at the layer top
      Rdir(:,k) = 0.5_dp*(sum_o(:) + dif_o(:))
      Tdir(:,k) = 0.5_dp*(sum_o(:) - dif_o(:))

    end do

    !! One-sided (Riccati) adding across the layers, with stream intensities

    ! Lower boundary - gas giant, no reflection
    RkN(:,nlev) = 0.0_dp
    RbkN(:,:,nlev) = 0.0_dp

    !! Direct beam transmission to level (running product of the layer transmissions)
    T(1) = 1.0_dp
    do k = 1, nlay
      T(k+1) = T(k)*dtr_z(k)
    end do

    !! Upward pass: reflection of everything below (RbkN) and beam source (RkN, per unit beam at the
    !! level), storing Minv_l = (I - R_l Rb_{l+1})^-1 for the downward pass
    do l = nlay, 1, -1
      lp1 = l + 1

      TT(:,:) = E(:,:) - matmul(RbkN(:,:,lp1),Rdiff(:,:,l))
      TT(:,:) = matinv2(TT(:,:))

      CC(:,:) = matmul(Tdiff(:,:,l),TT(:,:))

      DD(:) = RkN(:,lp1)*dtr(l) + matmul(RbkN(:,:,lp1),Tdir(:,l))
      RkN(:,l) = Rdir(:,l) + matmul(CC(:,:), DD(:))

      CC(:,:) = matmul(CC(:,:),RbkN(:,:,lp1))

      RbkN(:,:,l) = Rdiff(:,:,l) + matmul(CC(:,:),Tdiff(:,:,l))

      Minv(:,:,l) = matinv2(E(:,:) - matmul(Rdiff(:,:,l),RbkN(:,:,lp1)))

    end do

    !! Downward pass - upper boundary, no incident diffuse intensity
    D_I(:,1) = 0.0_dp
    do l = 1, nlay
      DD(:) = matmul(Tdiff(:,:,l),D_I(:,l)) + Tdir(:,l)*T(l) + matmul(Rdiff(:,:,l),RkN(:,l+1))*T(l+1)
      D_I(:,l+1) = matmul(Minv(:,:,l),DD(:))
    end do
    do k = 1, nlev
      U_I(:,k) = matmul(RbkN(:,:,k),D_I(:,k)) + RkN(:,k)*T(k)
    end do

    !! Fluxes - diffuse from the double-Gauss quadrature plus the direct beam, multiplied by the incident flux
    flx_down(:) = 0.0_dp
    flx_up(:) = 0.0_dp
    do i = 1, nmu
      flx_down(:) = flx_down(:) + D_I(i,:)*wa(i)*uarr(i)
      flx_up(:) = flx_up(:) + U_I(i,:)*wa(i)*uarr(i)
    end do
    flx_down(:) = (twopi*flx_down(:) + mu_z*T(:))*F0_in
    flx_up(:) = twopi*flx_up(:)*F0_in

  end subroutine sw_Liou88_method

  !! Legendre polynomials P_0..P_3
  pure function Pleg(l, x) result(p)
    implicit none
    integer, intent(in) :: l
    real(dp), intent(in) :: x
    real(dp) :: p
    select case(l)
    case(0)
      p = 1.0_dp
    case(1)
      p = x
    case(2)
      p = 0.5_dp*(3.0_dp*x**2 - 1.0_dp)
    case default
      p = 0.5_dp*(5.0_dp*x**3 - 3.0_dp*x)
    end select
  end function Pleg

  !! Eigenvector of a 2x2 matrix for eigenvalue lam, taking the better conditioned of the two
  !! row forms. If C is (near) a multiple of the identity any vector works, use the unit vector m.
  pure subroutine eigvec2(C, lam, m, v)
    implicit none
    real(dp), dimension(2,2), intent(in) :: C
    real(dp), intent(in) :: lam
    integer, intent(in) :: m
    real(dp), dimension(2), intent(out) :: v
    real(dp), dimension(2) :: v1, v2
    real(dp) :: n1, n2

    v1(1) = C(1,2) ; v1(2) = lam - C(1,1)
    v2(1) = lam - C(2,2) ; v2(2) = C(2,1)
    n1 = sqrt(v1(1)**2 + v1(2)**2)
    n2 = sqrt(v2(1)**2 + v2(2)**2)

    if (max(n1,n2) <= 1.0e-12_dp*max(abs(lam),1.0_dp)) then
      v(:) = 0.0_dp
      v(m) = 1.0_dp
    else if (n1 >= n2) then
      v(:) = v1(:)/n1
    else
      v(:) = v2(:)/n2
    end if

  end subroutine eigvec2

  pure function matinv2(A) result(B)
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
    real(dp), intent(in) :: A(2,2)   !! Matrix
    real(dp)             :: B(2,2)   !! Inverse matrix
    real(dp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1.0_dp/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = detinv * A(1,1)

  end function matinv2

end module sw_Liou88_mod
