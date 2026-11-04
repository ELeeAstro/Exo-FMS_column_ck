!!!
! Elspeth KH Lee - Aug 2023 : Initial version
!                - Nov 2026 : Bug fixes (upward flux sign, P coefficient, diffuse R matrix, upward sweep),
!                             delta-M (f = chi_4) as sw_SH4, TTHG option, single zenith angle, resonance guard
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!                - Nov 2026 : Optimisations - symmetric layer block reduction (2x2 inverses only), layer-last
!                             array layout, direct beam as a running product of the layer
!                             transmissions, reciprocals, contiguous per g-ordinance inputs
!
! sw: alpha-4SDA - 4 stream analytic spherical harmonic method with doubling adding method (Zhang & Li 2013)
!     Gives the same answer as sw_SH4 (to machine precision), but with 2x2 adding instead of a banded matrix
! Pros: Fast and accurate method with good multiple scattering approximation, no large matrix inversions
! Cons: Four-stream approximation
!!!

module sw_SDA_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Single scattering albedo cap - exactly conservative scattering (w = 1, eigenvalue k = 0) is singular.
  !! 1 - 1e-12 absorbs < 1e-8 of the incident flux even in a tau ~ 6000 column
  real(dp), parameter :: w_max = 1.0_dp - 1.0e-12_dp

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: fourpi = 4.0_dp * pi
  real(dp), parameter :: r_fourpi = 1.0_dp / fourpi

  !! Use Two-Term HG function (TTHG) for sph harmonic method, parameters are the PICASO defaults
  !! forward fraction = frac_a + frac_b*g_back**frac_c, g_forward = c_forward*g, g_back = c_back*g
  logical, parameter :: TTHG = .False.
  real(dp), parameter :: tthg_forward = 1.0_dp, tthg_back = -0.5_dp
  real(dp), parameter :: tthg_frac_a = 1.0_dp, tthg_frac_b = -1.0_dp
  integer, parameter :: tthg_frac_c = 2

  !! How close to the k*mu0 = 1 resonance before mu0 is nudged
  real(dp), parameter :: res_eps = 1.0e-6_dp

  private :: sw_fs_SDA, matinv2, matinv4, ludcmp, lubksb, inv_LU
  public :: sw_SDA

contains

  subroutine sw_SDA(nlay, nlev, nb, ng, gw, tau_e, mu_z, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
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
    integer :: b, g
    real(dp), dimension(nlev) :: sw_up_g, sw_down_g
    real(dp), dimension(nlev) :: tau_g
    real(dp), dimension(nlay) :: ssa_g, gg_g

    sw_up(:) = 0.0_dp
    sw_down(:) = 0.0_dp

    !! Shortwave flux calculation - sum over each band and g-ordinance
    if (mu_z > 0.0_dp) then
      do b = 1, nb
        if (Finc(b) <= 0.0_dp) cycle
        do g = 1, ng
          tau_g(:) = tau_e(g,b,:)
          ssa_g(:) = ssa(g,b,:)
          gg_g(:) = gg(g,b,:)
          call sw_fs_SDA(nlay, nlev, Finc(b), mu_z, tau_g, ssa_g, gg_g, sw_down_g, sw_up_g)
          sw_up(:) = sw_up(:) + sw_up_g(:) * gw(g)
          sw_down(:) = sw_down(:) + sw_down_g(:) * gw(g)
        end do
      end do
    end if

    !! Net sw flux
    sw_net(:) = sw_up(:) - sw_down(:)

    !! Absorbed Stellar Radiation (ASR)
    asr = sw_down(1) - sw_up(1)

  end subroutine sw_SDA

  subroutine sw_fs_SDA(nlay, nlev, F0_in, mu_in, tau_in, w_in, g_in, flx_down, flx_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0_in, mu_in
    real(dp), dimension(nlev), intent(in) :: tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    !! Work variables
    integer :: k
    real(dp), dimension(nlay) :: w0, dtau, g_f, g_b, f_f
    real(dp), dimension(nlev) :: T
    real(dp) :: f0
    real(dp) :: om0, om1, om2, om3
    real(dp) :: a0, a1, a2, a3, b0, b1, b2, b3
    real(dp) :: e1, e2
    real(dp) :: beta, gam, k1, k2, R1, R2, P1, P2, Q1, Q2
    real(dp) :: eta0, eta1, eta2, eta3, del0, del1, del2, del3, delta
    real(dp) :: z1p, z1m, z2p, z2m
    real(dp) :: phi1p, phi1m, phi2p, phi2m
    real(dp) :: Cphi1p, Cphi1m, Cphi2p, Cphi2m

    real(dp), dimension(2) :: h1t, h1b, hs, hd, Xs, Yd
    real(dp), dimension(2,2) :: Am, Ap, S1, D1, S2, D2, X, Y

    real(dp), dimension(2,nlay) :: Rdir, Tdir
    real(dp), dimension(2,2,nlay) :: Rdiff, Tdiff

    real(dp), dimension(2,nlev) :: RkN, D
    real(dp), dimension(2,2,nlev) :: RbkN
    real(dp), dimension(2,2,nlay) :: Minv

    real(dp), dimension(2,2) :: E, TT, CC
    real(dp), dimension(2) ::  DD

    integer, parameter :: nstr = 4
    real(dp), dimension(nlay) :: fc
    real(dp), dimension(nstr,nlay) :: chi
    real(dp), parameter :: eps_20 = 1.0e-20_dp

    integer :: l, lp1

    real(dp), dimension(nlay) :: dtr, dtr_z

    real(dp) :: mu_z, mu_l, rfc, rk1, rk2, ra3, rdelta, sq


    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w_in(:) <= 1.0e-12_dp)) then

      flx_down(:) = F0_in*mu_in * exp(-(tau_in(:) - tau_in(1))/mu_in)

      flx_up(:) = 0.0_dp

      return

    end if

    !! Calculate dtau in each layer
    dtau(:) = tau_in(2:) - tau_in(1:nlay)

    !! Phase function Legendre moments chi_l (HG: chi_l = g**l)
    if (TTHG .eqv. .False.) then
      chi(1,:) = g_in(:)
      do l = 2, nstr
        chi(l,:) = chi(l-1,:)*g_in(:)
      end do
    else
      g_f(:) = tthg_forward*g_in(:)
      g_b(:) = tthg_back*g_in(:)
      f_f(:) = tthg_frac_a + tthg_frac_b*g_b(:)**tthg_frac_c
      do l = 1, nstr
        chi(l,:) = f_f(:)*g_f(:)**l + (1.0_dp - f_f(:))*g_b(:)**l
      end do
    end if

    !! Delta-M scaling (f = chi_4, g**4 for HG as sw_SH4 and PICASO)
    fc(:) = chi(nstr,:)
    w0(:) = w_in(:)*((1.0_dp - fc(:))/(1.0_dp - fc(:)*w_in(:)))
    w0(:) = min(w0(:), w_max)
    dtau(:) = (1.0_dp - w_in(:)*fc(:))*dtau(:)

    !! Single zenith angle for the whole column (value at the lower boundary)
    mu_z = mu_in

    !! Start SDA calculation

    !! First find the Reflection and Transmission coefficents (direct and diffuse) for each layer
    do k = 1, nlay

      !! Omega Legendre polynomial coefficents - scale with delta-M (isotropic for g = 0)
      rfc = 1.0_dp/(1.0_dp - fc(k))
      om0 = 1.0_dp
      om1 = 3.0_dp * (chi(1,k) - fc(k))*rfc
      om2 = 5.0_dp * (chi(2,k) - fc(k))*rfc
      om3 = 7.0_dp * (chi(3,k) - fc(k))*rfc

      !! Find the a coefficents
      a0 =  1.0_dp - w0(k)*om0 + eps_20
      a1 =  3.0_dp - w0(k)*om1 + eps_20
      a2 =  5.0_dp - w0(k)*om2 + eps_20
      a3 =  7.0_dp - w0(k)*om3 + eps_20

      !! Find beta and gamma
      beta = a0*a1 + (4.0_dp/9.0_dp)*a0*a3 + (1.0_dp/9.0_dp)*a2*a3
      gam = (1.0_dp/9.0_dp)*a0*a1*a2*a3

      !! Find k values - lambda in Rooney
      sq = sqrt(beta**2 - 4.0_dp*gam)
      k1 = sqrt(0.5_dp*(beta + sq))
      k2 = sqrt(0.5_dp*(beta - sq))
      rk1 = 1.0_dp/k1
      rk2 = 1.0_dp/k2

      !! Zenith angle for this layer - nudge mu0 off the k*mu0 = 1 resonance (0/0 in the particular solution)
      mu_l = mu_z
      if (abs(1.0_dp - k1*mu_l) < res_eps) mu_l = mu_l + 2.0_dp*res_eps*rk1
      if (abs(1.0_dp - k2*mu_l) < res_eps) mu_l = mu_l + 2.0_dp*res_eps*rk2

      !! Inverse zenith angle
      f0 = 1.0_dp/mu_l

      !! Layer transmission (and at the un-nudged mu0 for the direct beam)
      dtr(k) = exp(-dtau(k)*f0)
      if (mu_l == mu_z) then
        dtr_z(k) = dtr(k)
      else
        dtr_z(k) = exp(-dtau(k)/mu_z)
      end if

      !! Find the b coefficents - normalise Finc to 1 here
      b0 = w0(k)*om0 * r_fourpi
      b1 = w0(k)*om1 * -(mu_l) * r_fourpi
      b2 = 0.5_dp * w0(k)*om2 * (3.0_dp * mu_l**2 - 1.0_dp) * r_fourpi
      b3 = 0.5_dp * w0(k)*om3 * (5.0_dp * -mu_l**3 - 3.0_dp*-(mu_l)) * r_fourpi

      !! Find e values
      e1 = exp(-k1*dtau(k))
      e2 = exp(-k2*dtau(k))      
      
      !! Find R, P and Q coefficents 
      !! NOTE: Zhang et al. (2013) has the wrong coefficent definitions
      !! Rooney et al. (2023) has the correct definitions and order in the matrix
      !! So we use the Rooney definitions, but keep the Zhang notation
      ra3 = 1.0_dp/a3
      Q1 = -1.5_dp * (a0*a1*rk1 - k1)*ra3
      Q2 = -1.5_dp * (a0*a1*rk2 - k2)*ra3
      R1 = -a0*rk1
      R2 = -a0*rk2
      P1 = 0.5_dp * (a0*a1*rk1**2 - 1.0_dp)
      P2 = 0.5_dp * (a0*a1*rk2**2 - 1.0_dp)

      !! Find the delta values
      delta = 9.0_dp*(f0**4 - beta*f0**2 + gam)
      del0 = (a1*b0 - b1*f0)*(a2*a3 - 9.0_dp*f0**2) + 2.0_dp*f0**2*(a3*b2 - 2.0_dp*a3*b0 - 3.0_dp*b3*f0)
      del1 = (a0*b1 - b0*f0)*(a2*a3 - 9.0_dp*f0**2) - 2.0_dp*a0*f0*(a3*b2 - 3.0_dp*b3*f0)
      del2 = (a3*b2 - 3.0_dp*b3*f0)*(a0*a1 - f0**2) - 2.0_dp*a3*f0*(a0*b1 - b0*f0)
      del3 = (a2*b3 - 3.0_dp*b2*f0)*(a0*a1 - f0**2) + f0**2*(6.0_dp*a0*b1 - 4.0_dp*a0*b3 - 6.0_dp*b0*f0)

      !! Find the eta values
      rdelta = 1.0_dp/delta
      eta0 = del0*rdelta
      eta1 = del1*rdelta
      eta2 = 0.625_dp * del2*rdelta
      eta3 = del3*rdelta

      !! Find the phi values
      phi1p = twopi*(0.5_dp + R1 + 5.0_dp/8.0_dp*P1)
      phi1m = twopi*(0.5_dp - R1 + 5.0_dp/8.0_dp*P1)
      phi2p = twopi*(0.5_dp + R2 + 5.0_dp/8.0_dp*P2)
      phi2m = twopi*(0.5_dp - R2 + 5.0_dp/8.0_dp*P2)

      !! Find the Phi values
      Cphi1p = twopi*(-1.0_dp/8.0_dp + 5.0_dp/8.0_dp*P1 + Q1)
      Cphi1m = twopi*(-1.0_dp/8.0_dp + 5.0_dp/8.0_dp*P1 - Q1)
      Cphi2p = twopi*(-1.0_dp/8.0_dp + 5.0_dp/8.0_dp*P2 + Q2)
      Cphi2m = twopi*(-1.0_dp/8.0_dp + 5.0_dp/8.0_dp*P2 - Q2)

      !! Find the Z values
      z1p = twopi*(0.5_dp*eta0 + eta1 + eta2)
      z1m = twopi*(0.5_dp*eta0 - eta1 + eta2)
      z2p = twopi*(-1.0_dp/8.0_dp*eta0 + eta2 + eta3)
      z2m = twopi*(-1.0_dp/8.0_dp*eta0 + eta2 - eta3)

      !! Boundary value problem. With the columns ordered (mode 1 and 2 decaying from the top | from
      !! the bottom), the 4x4 matrices of Zhang & Li (2013) / Rooney et al. (2023) have the layer's
      !! up-down symmetry, A1 = [[Am, Ap E], [Ap E, Am]] and A2 = [[Am E, Ap], [Ap, Am E]], E = diag(e1, e2),
      !! rows (phi, Phi) at the two boundaries. Both are block diagonal in the sum / difference basis, so
      !! A2 A1^-1 = 1/2 [[X + Y, X - Y], [X - Y, X + Y]],  X = S2 S1^-1,
      !! Y = D2 D1^-1, S1/D1 = Am +/- Ap E, S2/D2 = Am E +/- Ap - two 2x2 inverses instead of the 4x4
      !! adjugate inverse.
      Am(1,1) = phi1m ; Am(1,2) = phi2m ; Am(2,1) = Cphi1m ; Am(2,2) = Cphi2m
      Ap(1,1) = phi1p ; Ap(1,2) = phi2p ; Ap(2,1) = Cphi1p ; Ap(2,2) = Cphi2p
      S1(:,1) = Am(:,1) + Ap(:,1)*e1 ; S1(:,2) = Am(:,2) + Ap(:,2)*e2
      D1(:,1) = Am(:,1) - Ap(:,1)*e1 ; D1(:,2) = Am(:,2) - Ap(:,2)*e2
      S2(:,1) = Am(:,1)*e1 + Ap(:,1) ; S2(:,2) = Am(:,2)*e2 + Ap(:,2)
      D2(:,1) = Am(:,1)*e1 - Ap(:,1) ; D2(:,2) = Am(:,2)*e2 - Ap(:,2)
      X(:,:) = matmul(S2(:,:), matinv2(S1(:,:)))
      Y(:,:) = matmul(D2(:,:), matinv2(D1(:,:)))

      !! Diffuse reflection and transmission (no normalisation): Tdiff = (X + Y)/2, Rdiff = (X - Y)/2
      Tdiff(:,:,k) = 0.5_dp*(X(:,:) + Y(:,:))
      Rdiff(:,:,k) = 0.5_dp*(X(:,:) - Y(:,:))

      !! Direct beam: F = A2 A1^-1 H1 + H2, with H1 = (-z1m, -z2m | -z1p T, -z2p T),
      !! H2 = (z1m T, z2m T | z1p, z2p) - normalised by the beam flux (= 1 here)
      h1t(1) = -z1m ; h1t(2) = -z2m
      h1b(1) = -z1p*dtr(k) ; h1b(2) = -z2p*dtr(k)
      hs(:) = h1t(:) + h1b(:)
      hd(:) = h1t(:) - h1b(:)
      Xs(:) = matmul(X(:,:), hs(:))
      Yd(:) = matmul(Y(:,:), hd(:))
      Tdir(1,k) = (0.5_dp*(Xs(1) + Yd(1)) + z1m*dtr(k))*f0
      Tdir(2,k) = (0.5_dp*(Xs(2) + Yd(2)) + z2m*dtr(k))*f0
      Rdir(1,k) = (0.5_dp*(Xs(1) - Yd(1)) + z1p)*f0
      Rdir(2,k) = (0.5_dp*(Xs(2) - Yd(2)) + z2p)*f0

    end do

    !! We now have the transmission and reflection coefficents for both the direct and diffuse components
    !! Now we perform the doubling-adding method accros multiple layers

    !! One-sided (Riccati) adding: an upward pass for the reflection of everything below each level
    !! (RbkN) and the beam source (RkN), storing Minv_l = (I - R_l Rb_{l+1})^-1, then one downward sweep
    !! D_{l+1} = Minv_l [Tdiff_l D_l + Tdir_l T_l + Rdiff_l RkN_{l+1} T_{l+1}],  U_k = RbkN_k D_k + RkN_k T_k
    !! (identical to the two-sided adding, ~1/3 less work)

    ! Lower boundary - gas giant, no reflection
    RkN(:,nlev) = 0.0_dp
    RbkN(:,:,nlev) = 0.0_dp

    !! Direct beam transmission to level - running product of the layer transmissions (a new
    !! exponential only where mu0 was nudged off a resonance in that layer)
    T(1) = 1.0_dp
    do k = 1, nlay
      T(k+1) = T(k)*dtr_z(k)
    end do

    !! E indentity matrix
    E(1,1) = 1.0_dp ; E(1,2) = 0.0_dp ; E(2,1) = 0.0_dp ; E(2,2) = 1.0_dp

    do l = nlay, 1, -1
      lp1 = l + 1

      TT(:,:) = matinv2(E(:,:) - matmul(RbkN(:,:,lp1),Rdiff(:,:,l)))

      CC(:,:) = matmul(Tdiff(:,:,l),TT(:,:))

      DD(:) = RkN(:,lp1)*dtr(l) + matmul(RbkN(:,:,lp1),Tdir(:,l))
      RkN(:,l) = Rdir(:,l) + matmul(CC(:,:), DD(:))

      CC(:,:) = matmul(CC(:,:),RbkN(:,:,lp1))

      RbkN(:,:,l) = Rdiff(:,:,l) + matmul(CC(:,:),Tdiff(:,:,l))

      Minv(:,:,l) = matinv2(E(:,:) - matmul(Rdiff(:,:,l),RbkN(:,:,lp1)))

    end do

    D(:,1) = 0.0_dp
    do l = 1, nlay
      DD(:) = matmul(Tdiff(:,:,l),D(:,l)) + Tdir(:,l)*T(l) + matmul(Rdiff(:,:,l),RkN(:,l+1))*T(l+1)
      D(:,l+1) = matmul(Minv(:,:,l),DD(:))
    end do

    !! Down and up fluxes are multiplied by the incident flux (U_k = RbkN_k D_k + RkN_k T_k, first moment)
    do k = 1, nlev
      flx_down(k) = (D(1,k) + T(k))*mu_z*F0_in
      flx_up(k) = (RbkN(1,1,k)*D(1,k) + RbkN(1,2,k)*D(2,k) + RkN(1,k)*T(k))*mu_z*F0_in
    end do

  end subroutine sw_fs_SDA

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

  pure function matinv4(A) result(B)
    !! Performs a direct calculation of the inverse of a 4×4 matrix.
    real(dp), intent(in) :: A(4,4)   !! Matrix
    real(dp)             :: B(4,4)   !! Inverse matrix
    real(dp)             :: detinv, s0, s1, s2, s3, s4, s5, c5, c4, c3, c2, c1, c0

    s0 = A(1,1) * A(2,2) - A(2,1) * A(1,2)
    s1 = A(1,1) * A(2,3) - A(2,1) * A(1,3)
    s2 = A(1,1) * A(2,4) - A(2,1) * A(1,4)
    s3 = A(1,2) * A(2,3) - A(2,2) * A(1,3)
    s4 = A(1,2) * A(2,4) - A(2,2) * A(1,4)
    s5 = A(1,3) * A(2,4) - A(2,3) * A(1,4)

    c5 = A(3,3) * A(4,4) - A(4,3) * A(3,4)
    c4 = A(3,2) * A(4,4) - A(4,2) * A(3,4)
    c3 = A(3,2) * A(4,3) - A(4,2) * A(3,3)
    c2 = A(3,1) * A(4,4) - A(4,1) * A(3,4)
    c1 = A(3,1) * A(4,3) - A(4,1) * A(3,3)
    c0 = A(3,1) * A(4,2) - A(4,1) * A(3,2)

    detinv = 1.0_dp / (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0)

    B(1,1) = ( A(2,2) * c5 - A(2,3) * c4 + A(2,4) * c3) * detinv
    B(1,2) = (-A(1,2) * c5 + A(1,3) * c4 - A(1,4) * c3) * detinv
    B(1,3) = ( A(4,2) * s5 - A(4,3) * s4 + A(4,4) * s3) * detinv
    B(1,4) = (-A(3,2) * s5 + A(3,3) * s4 - A(3,4) * s3) * detinv

    B(2,1) = (-A(2,1) * c5 + A(2,3) * c2 - A(2,4) * c1) * detinv
    B(2,2) = ( A(1,1) * c5 - A(1,3) * c2 + A(1,4) * c1) * detinv
    B(2,3) = (-A(4,1) * s5 + A(4,3) * s2 - A(4,4) * s1) * detinv
    B(2,4) = ( A(3,1) * s5 - A(3,3) * s2 + A(3,4) * s1) * detinv

    B(3,1) = ( A(2,1) * c4 - A(2,2) * c2 + A(2,4) * c0) * detinv
    B(3,2) = (-A(1,1) * c4 + A(1,2) * c2 - A(1,4) * c0) * detinv
    B(3,3) = ( A(4,1) * s4 - A(4,2) * s2 + A(4,4) * s0) * detinv
    B(3,4) = (-A(3,1) * s4 + A(3,2) * s2 - A(3,4) * s0) * detinv

    B(4,1) = (-A(2,1) * c3 + A(2,2) * c1 - A(2,3) * c0) * detinv
    B(4,2) = ( A(1,1) * c3 - A(1,2) * c1 + A(1,3) * c0) * detinv
    B(4,3) = (-A(4,1) * s3 + A(4,2) * s1 - A(4,3) * s0) * detinv
    B(4,4) = ( A(3,1) * s3 - A(3,2) * s1 + A(3,3) * s0) * detinv

  end function matinv4

  subroutine ludcmp(A,n,np,indx,D)
    implicit none

    integer, intent(in) :: n, np
    real(dp), dimension(np,np), intent(inout) :: A

    integer, dimension(n), intent(out) :: indx
    real(dp), intent(out) :: D

    integer, parameter :: nmax = 100
    real(dp), parameter :: tiny = 1.0e-20_dp
    real(dp), dimension(nmax) :: vv

    integer :: i, j, k, imax
    real(dp) :: aamax, dum, sum

    D = 1.0_dp

    do i = 1, n
      aamax = 0.0_dp
      do j = 1, n
        if (abs(A(i,j)) > aamax) then
          aamax = abs(A(i,j))
        end if
      end do
      if (aamax == 0.0_dp) then
        print*, 'singualr matrix in LU decomp!'
        stop
      end if
      vv(i) = 1.0_dp/aamax
    end do

  
    do j = 1, n
      do i = 1, j-1
        sum = A(i,j)
        do k = 1, i-1
          sum = sum  - A(i,k)*A(k,j)
        end do
        A(i,j) = sum
      end do
      aamax = 0.0_dp
      do i = j, n
        sum = A(i,j)
        do k = 1, j-1
          sum = sum  - A(i,k)*A(k,j)
        end do
        A(i,j) = sum
        dum = vv(i)*abs(sum)
        if (dum >= aamax) then
          imax = i
          aamax = dum
        end if
      end do
      if (j /= imax) then
        do k = 1, n
          dum = A(imax,k)
          A(imax,k) = A(j,k)
          A(j,k) = dum
        end do 
        D = -D
        vv(imax) = vv(j)
      end if
      indx(j) = imax
      if (A(j,j) <= tiny) then
        A(j,j) = tiny
      end if
      if (j /= n) then
        dum = 1.0_dp/A(j,j)
        do i = j+1, n
          A(i,j) = A(i,j)*dum
        end do
      end if
    end do

  end subroutine ludcmp

  subroutine lubksb(A, n, np, indx, B)
    implicit none

    integer, intent(in) :: n, np
    integer, dimension(n), intent(in) :: indx
    real(dp), dimension(np,np), intent(in) :: A

    real(dp), dimension(n), intent(out) :: B

    integer :: i, j, ii, ll
    real(dp) :: sum

    ii = 0

    do i = 1, n
      ll = indx(i)
      sum = B(ll)
      b(ll) = b(i)
      if (ii /= 0) then
        do j = ii,i-1
          sum = sum - A(i,j)*B(j)
        end do
      else if (sum /= 0.0_dp) then
        ii = i
      end if
      B(i) = sum
    end do

    do i = n, 1, -1
      sum = B(i)
      if (i < n) then
        do j = i+1, n
          sum = sum - A(i,j)*B(j)
        end do
      end if
      B(i) = sum/A(i,i)
    end do
    
  end subroutine lubksb

  subroutine inv_LU(A,n,np,Y)
    implicit none

    integer, intent(in) :: n, np
    real(dp), dimension(np,np), intent(inout) :: A

    integer, dimension(n) :: indx
    real(dp), dimension(np,np), intent(out) :: Y

    real(dp) :: D
    integer :: i, j

    do i = 1, n
      do j = 1, n
        Y(i,j) = 0.0_dp
      end do
      Y(i,i) = 1.0_dp
    end do

    call ludcmp(A,n,np,indx,D)

    do j = 1, n
      call lubksb(A,n,np,indx,Y(1,j))
    end do

  end subroutine inv_LU

end module sw_SDA_mod
