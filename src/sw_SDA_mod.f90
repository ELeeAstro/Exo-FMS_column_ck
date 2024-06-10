module sw_SDA_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: fourpi = 4.0_dp * pi

  !! Use Two-Term HG function for sph harmonic method
  logical, parameter :: TTHG = .False.

  private :: sw_fs_SDA, matinv2, matinv4, ludcmp, lubksb, inv_LU
  public :: sw_SDA

contains

  subroutine sw_SDA(nlay, nlev, nb, ng, gw, tau_e, mu_z, Finc, ssa, gg, a_surf, sw_up, sw_down, sw_net, asr)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev, nb, ng
    real(dp), dimension(ng) :: gw
    real(dp), dimension(ng,nb,nlev), intent(in) :: tau_e
    real(dp), dimension(ng,nb,nlay), intent(in) :: ssa, gg
    real(dp), dimension(nlev) :: mu_z
    real(dp), dimension(nb) :: Finc, a_surf

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: sw_up, sw_down, sw_net
    real(dp), intent(out) :: asr

    !! Work variables
    integer :: b, g
    real(dp), dimension(nb,nlev) :: sw_down_b, sw_up_b
    real(dp), dimension(ng,nlev) :: sw_down_g, sw_up_g


    !! Shortwave flux calculation
    if (mu_z(nlev) > 0.0_dp) then
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
      do b = 1, nb
        sw_down_b(b,:) = 0.0_dp
        sw_up_b(b,:) = 0.0_dp
        if (Finc(b) <= 0.0_dp) then
          cycle
        end if
        do g = 1, ng
          call sw_fs_SDA(nlay, nlev, Finc(b), mu_z(:), tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), a_surf(b), sw_down_g(g,:), sw_up_g(g,:))
          sw_down_b(b,:) = sw_down_b(b,:) + sw_down_g(g,:) * gw(g)
          sw_up_b(b,:) = sw_up_b(b,:) + sw_up_g(g,:) * gw(g)
        end do
        sw_down(:) = sw_down(:) + sw_down_b(b,:)
        sw_up(:) = sw_up(:) + sw_up_b(b,:)
      end do
    else
      sw_down(:) = 0.0_dp
      sw_up(:) = 0.0_dp
    end if

    !! Net sw flux
    sw_net(:) = sw_up(:) - sw_down(:)

    !! Absorbed Stellar Radiation (ASR)
    asr = sw_down(1) - sw_up(1)

  end subroutine sw_SDA

  subroutine sw_fs_SDA(nlay, nlev, F0_in, mu_in, tau_in, w_in, g_in, w_surf_in, flx_down, flx_up)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: F0_in, w_surf_in
    real(dp), dimension(nlev), intent(in) :: tau_in, mu_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_down, flx_up

    !! Work variables
    integer :: k
    real(dp), dimension(nlay) :: w0, dtau, hg
    real(dp), dimension(nlev) :: tau, T
    real(dp) :: f0
    real(dp) :: om0, om1, om2, om3
    real(dp) :: a0, a1, a2, a3, b0, b1, b2, b3
    real(dp) :: e1, e2
    real(dp) :: beta, gam, k1, k2, R1, R2, P1, P2, Q1, Q2
    real(dp) :: eta0, eta1, eta2, eta3, del0, del1, del2, del3, delta
    real(dp) :: z1p, z1m, z2p, z2m
    real(dp) :: phi1p, phi1m, phi2p, phi2m
    real(dp) :: Cphi1p, Cphi1m, Cphi2p, Cphi2m

    real(dp), dimension(4) :: H1, H2, H3, H4
    real(dp), dimension(4) :: AA12H1
    real(dp), dimension(4,4) :: AA1, AA2, AA1_i, AA12
    real(dp), dimension(4) :: Fdir, Fdiffa, Fdiffb

    real(dp), dimension(nlay,2) :: Rdir, Tdir
    real(dp), dimension(nlay,2,2) :: Rdiff, Tdiff 

    real(dp), dimension(nlev,2) :: T1k, RkN, U, D
    real(dp), dimension(nlev,2,2) :: Rst1k, RbkN

    real(dp), dimension(2,2) :: E, TT, CC
    real(dp), dimension(2) ::  DD

    real(dp), dimension(nlay) :: fc, sigma_sq, pmom2, c
    integer, parameter :: nstr = 4
    real(dp), parameter :: eps_20 = 1.0e-20_dp

    integer :: l, km1, lp1

    real(dp), dimension(nlay) :: dtr

    real(dp) :: hg2, alp

    real(dp), dimension(nlev) :: cum_trans

    !! If zero albedo across all atmospheric layers then return direct beam only
    if (all(w_in(:) <= 1.0e-12_dp)) then

      ! Check mu_z for spherical corrections
      if (mu_in(nlev) == mu_in(1)) then
        ! No zenith correction, use regular method
        flx_down(:) = F0_in* mu_in(nlev) * exp(-tau_in(:)/mu_in(nlev))
      else
        ! Zenith angle correction, use cumulative transmission function
        cum_trans(1) = tau_in(1)/mu_in(1)
        do k = 1, nlev-1
          cum_trans(k+1) = cum_trans(k) + (tau_in(k+1) - tau_in(k))/mu_in(k+1)
        end do
        do k = 1, nlev
          flx_down(k) = F0_in * mu_in(nlev) * exp(-cum_trans(k))
        end do
      end if

      flx_down(nlev) = flx_down(nlev) * (1.0_dp - w_surf_in) ! The surface flux for surface heating is the amount of flux absorbed by surface
      flx_up(:) = 0.0_dp ! We assume no upward flux here even if surface albedo

      return

    end if

    !! Calculate dtau in each layer
    dtau(:) = tau_in(2:) - tau_in(1:nlay)

    !! Delta-M+ scaling (Following DISORT: Lin et al. 2018)
    !! Assume HG phase function for scaling

    where (g_in(:) >= 1.0e-6_dp)
      fc(:) = g_in(:)**(nstr)
      pmom2(:) = g_in(:)**(nstr+1)
      sigma_sq(:) = real((nstr+1)**2 - nstr**2,dp) / &
      & ( log(fc(:)**2/pmom2(:)**2) )
      c(:) = exp(real(nstr**2,dp)/(2.0_dp*sigma_sq(:)))
      fc(:) = c(:)*fc(:)

      w0(:) = w_in(:)*((1.0_dp - fc(:))/(1.0_dp - fc(:)*w_in(:)))
      dtau(:) = (1.0_dp - w_in(:)*fc(:))*dtau(:)

    elsewhere
      w0(:) = w_in(:)
      fc(:) = 0.0_dp
    end where

    hg(:) = g_in(:)

    !! Reform edge optical depths
    tau(1) = tau_in(1)
    do k = 1, nlay
      tau(k+1) = tau(k) + dtau(k)
    end do

    !! Start SDA calculation

    !! First find the Reflection and Transmission coefficents (direct and diffuse) for each layer
    do k = 1, nlay

      !! Layer transmission
      dtr(k) = exp(-dtau(k)/mu_in(k))

      !! Inverse zenith angle
      f0 = 1.0_dp/mu_in(k)

      !! Omega Legendre polynomial coefficents - scale with delta-M+
      if (hg(k) >= 1.0e-6_dp) then
        if (TTHG .eqv. .False.) then
          ! Use HG phase function
          om0 = 1.0_dp
          om1 = 3.0_dp * (hg(k) - fc(k))/(1.0_dp - fc(k))
          om2 = 5.0_dp * (hg(k)**2 - fc(k))/(1.0_dp - fc(k))
          om3 = 7.0_dp * (hg(k)**3 - fc(k))/(1.0_dp - fc(k))
        else
          ! Use TTHG phase function with default parameters
          hg2 = hg(k)/2.0_dp
          alp = 1.0_dp - hg2**2
          om0 = 1.0_dp
          om1 = 3.0_dp * ((alp*hg(k) + (1.0_dp - alp)*hg2) - fc(k))/(1.0_dp - fc(k))
          om2 = 5.0_dp * ((alp*hg(k)**2 + (1.0_dp - alp)*hg2**2) - fc(k))/(1.0_dp - fc(k))
          om3 = 7.0_dp * ((alp*hg(k)**3 + (1.0_dp - alp)*hg2**3) - fc(k))/(1.0_dp - fc(k))
        end if
      else
        ! Use Rayleigh scattering phase function for isotropic scattering
        om0 = 1.0_dp
        om1 = 0.0_dp
        om2 = 0.5_dp
        om3 = 0.0_dp
      end if

      !! Find the a coefficents
      a0 =  1.0_dp - w0(k)*om0 + eps_20
      a1 =  3.0_dp - w0(k)*om1 + eps_20
      a2 =  5.0_dp - w0(k)*om2 + eps_20
      a3 =  7.0_dp - w0(k)*om3 + eps_20

      !! Find the b coefficents - normalise Finc to 1 here
      b0 = w0(k)*om0 / fourpi
      b1 = w0(k)*om1 * -(mu_in(k)) / fourpi
      b2 = 0.5_dp * w0(k)*om2 * (3.0_dp * mu_in(k)**2 - 1.0_dp) / fourpi
      b3 = 0.5_dp * w0(k)*om3 * (5.0_dp * -mu_in(k)**3 - 3.0_dp*-(mu_in(k))) / fourpi

      !! Find beta and gamma
      beta = a0*a1 + (4.0_dp/9.0_dp)*a0*a3 + (1.0_dp/9.0_dp)*a2*a3
      gam = (1.0_dp/9.0_dp)*a0*a1*a2*a3

      !! Find k values - lambda in Rooney
      k1 = (beta + sqrt((beta**2 - 4.0_dp*gam)))/2.0_dp
      k2 = (beta - sqrt((beta**2 - 4.0_dp*gam)))/2.0_dp

      k1 = sqrt(k1)
      k2 = sqrt(k2)

      !! Find e values
      e1 = exp(-k1*dtau(k))
      e2 = exp(-k2*dtau(k))      
      
      !! Find R, P and Q coefficents 
      !! NOTE: Zhang et al. (2013) has the wrong coefficent definitions
      !! Rooney et al. (2023) has the correct definitions and order in the matrix
      !! So we use the Rooney definitions, but keep the Zhang notation
      Q1 = -3.0_dp/2.0_dp * (a0*a1/k1 - k1)/a3
      Q2 = -3.0_dp/2.0_dp * (a0*a1/k2 - k2)/a3
      R1 = -a0/k1
      R2 = -a0/k2
      P1 = 0.3125_dp * (a0*a1/k1**2 - 1.0_dp)
      P2 = 0.3125_dp * (a0*a1/k2**2 - 1.0_dp)

      !! Find the delta values
      delta = 9.0_dp*(f0**4 - beta*f0**2 + gam)
      del0 = (a1*b0 - b1*f0)*(a2*a3 - 9.0_dp*f0**2) + 2.0_dp*f0**2*(a3*b2 - 2.0_dp*a3*b0 - 3.0_dp*b3*f0)
      del1 = (a0*b1 - b0*f0)*(a2*a3 - 9.0_dp*f0**2) - 2.0_dp*a0*f0*(a3*b2 - 3.0_dp*b3*f0)
      del2 = (a3*b2 - 3.0_dp*b3*f0)*(a0*a1 - f0**2) - 2.0_dp*a3*f0*(a0*b1 - b0*f0)
      del3 = (a2*b3 - 3.0_dp*b2*f0)*(a0*a1 - f0**2) + f0**2*(6.0_dp*a0*b1 - 4.0_dp*a0*b3 - 6.0_dp*b0*f0)

      !! Find the eta values
      eta0 = del0/delta
      eta1 = del1/delta
      eta2 = 0.625_dp * del2/delta
      eta3 = del3/delta

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

      !! Find A1 matrix
      AA1(1,1) = phi1m ; AA1(1,2) = phi1p*e1 ; AA1(1,3) = phi2m ; AA1(1,4) = phi2p*e2
      AA1(2,1) = Cphi1m ; AA1(2,2) = Cphi1p*e1 ; AA1(2,3) = Cphi2m; AA1(2,4) = Cphi2p*e2
      AA1(3,1) = phi1p*e1 ; AA1(3,2) = phi1m ; AA1(3,3) = phi2p*e2 ; AA1(3,4) = phi2m    
      AA1(4,1) = Cphi1p*e1 ; AA1(4,2) = Cphi1m ; AA1(4,3) = Cphi2p*e2 ; AA1(4,4) = Cphi2m

      !! Find H1 vector
      H1(1) = -z1m ;  H1(2) = -z2m ;  H1(3) = -z1p * dtr(k) ; H1(4) = -z2p * dtr(k)

      !! Find A2 matrix
      AA2(1,1) = phi1m*e1 ; AA2(1,2) = phi1p ; AA2(1,3) = phi2m*e2 ; AA2(1,4) = phi2p
      AA2(2,1) = Cphi1m*e1 ; AA2(2,2) = Cphi1p ; AA2(2,3) = Cphi2m*e2; AA2(2,4) = Cphi2p
      AA2(3,1) = phi1p ; AA2(3,2) = phi1m*e1 ; AA2(3,3) = phi2p ; AA2(3,4) = phi2m*e2    
      AA2(4,1) = Cphi1p ; AA2(4,2) = Cphi1m*e1 ; AA2(4,3) = Cphi2p ; AA2(4,4) = Cphi2m*e2

      !! Find H2 vector
      H2(1) = z1m * dtr(k) ;  H2(2) = z2m * dtr(k) ;  H2(3) = z1p ; H2(4) = z2p

      !! Now we need to invert the A1 matrix - Zhang and Li (2013) use a adjugate matrix method 
      !! with some reduction in the matrix order or used symetrical term (not 100% sure what they did)
      !! We primarily use the same method but keep the 4x4 layout
      !! We have LU decomposition here as an alternative in case of numerical instability (and for testing)

      AA1_i(:,:) = matinv4(AA1(:,:)) ! Use matrix determinant and adjugate method (faster but can be numerically unstable)
      !call inv_LU(AA1,4,4,AA1_i)    ! Use LU decomposition (slower but probably more stable)

      !! Multiply the AA1_i and AA2 matrices
      AA12(:,:) = matmul(AA2(:,:),AA1_i(:,:))

      !! Multiply the AA12 and H2 matrix
      AA12H1(:) = matmul(AA12(:,:),H1(:))

      !! Direct flux component - now we have the flux array (F(1)(2) = neg flux at lower,F(3)(4) = pos flux at upper)
      Fdir(:) = AA12H1(:) + H2(:)

      !! Store the direct beam reflection and transmission for this layer - normalised by the beam flux (=1 here)
      Rdir(k,1) = Fdir(3)/(mu_in(k))
      Rdir(k,2) = Fdir(4)/(mu_in(k))

      Tdir(k,1) = Fdir(1)/(mu_in(k))
      Tdir(k,2) = Fdir(2)/(mu_in(k))

      !! Now find the diffusive flux component

      !! Vector H3 is
      H3(1) = 1.0_dp; H3(2) = 0.0_dp; H3(3) = 0.0_dp; H3(4) = 0.0_dp 

      !! Multiply the AA12 and H3 matrix to get `a' boundary fluxes at layer edges (levels)
      Fdiffa(:) = matmul(AA12(:,:),H3(:))

      !! Vector H4 is
      H4(1) = 0.0_dp; H4(2) = 1.0_dp; H4(3) = 0.0_dp; H4(4) = 0.0_dp 

      !! Multiply the AA12 and H4 matrix to get `b' boundary fluxes at layer edges (levels)
      Fdiffb(:) = matmul(AA12(:,:),H4(:))

      !! Store the diffuse reflection and transmission for this layer - no normalisation
      Rdiff(k,1,1) = Fdiffa(3)
      Rdiff(k,1,2) = Fdiffb(3)
      Rdiff(k,2,2) = Fdiffa(4)
      Rdiff(k,2,1) = Fdiffb(4)
  
      Tdiff(k,1,1) = Fdiffa(1)
      Tdiff(k,1,2) = Fdiffb(1)
      Tdiff(k,2,1) = Fdiffa(2)
      Tdiff(k,2,2) = Fdiffb(2)

    end do

    !! We now have the transmission and reflection coefficents for both the direct and diffuse components
    !! Now we perform the doubling-adding method accros multiple layers

    !! Do boundary conditons first

    ! Upper
    T1k(1,:) = 0.0_dp
    Rst1k(1,:,:) = 0.0_dp

    ! Lower
    RkN(nlev,1) = w_surf_in ; RkN(nlev,2) = -w_surf_in/4.0_dp 
    RbkN(nlev,1,1) = w_surf_in ; RbkN(nlev,1,2) = 0.0_dp
    RbkN(nlev,2,1) = -w_surf_in/4.0_dp ; RbkN(nlev,2,2) = 0.0_dp

    !! Direct beam transmission to level
    T(:) = exp(-tau(:)/mu_in(k))

    !! E indentity matrix
    E(1,1) = 1.0_dp ; E(1,2) = 0.0_dp ; E(2,1) = 0.0_dp ; E(2,2) = 1.0_dp

    do k = 2, nlev
      km1 = k - 1        

      TT(:,:) = E(:,:) - matmul(Rdiff(km1,:,:),Rst1k(km1,:,:))
      TT(:,:) = matinv2(TT(:,:))

      CC(:,:) = matmul(Tdiff(km1,:,:),Rst1k(km1,:,:))
      CC(:,:) = matmul(CC(:,:),TT(:,:))

      DD(:) = Rdir(km1,:)*T(km1) + matmul(Rdiff(km1,:,:),T1k(km1,:))

      T1k(k,:) = Tdir(km1,:)*T(km1) + matmul(Tdiff(km1,:,:), T1k(km1,:)) + matmul(CC(:,:),DD(:))

      Rst1k(k,:,:) = Rdiff(km1,:,:) + matmul(CC(:,:),Tdiff(km1,:,:))

    end do
  
    do l = nlay, 1, -1
      lp1 = l + 1

      TT(:,:) = E(:,:) - matmul(RbkN(lp1,:,:),Rdiff(l,:,:))
      TT(:,:) = matinv2(TT(:,:))

      CC(:,:) = matmul(Tdiff(l,:,:),TT(:,:))

      DD(:) = RkN(lp1,:)*dtr(l) + matmul(RbkN(lp1,:,:),Tdir(l,:))
      DD(:) = matmul(CC(:,:),DD(:))

      RkN(l,:) = Rdir(l,:) + matmul(CC(:,:), DD(:))

      CC(:,:) = matmul(CC(:,:),RbkN(lp1,:,:))

      RbkN(l,:,:) = Rdiff(l,:,:) + matmul(CC(:,:),Tdiff(l,:,:))

    end do

    do k = 1, nlev

      TT(:,:) = E(:,:) - matmul(RbkN(k,:,:),Rst1k(k,:,:))
      TT(:,:) = matinv2(TT(:,:))

      DD(:) = RkN(k,:)*T(k) + matmul(RbkN(k,:,:),T1k(k,:))

      U(k,:) = matmul(TT(:,:),DD(:))

      D(k,:) = T1k(k,:) + matmul(Rst1k(k,:,:),U(k,:))

    end do

    !! Down and up fluxes are multiplied by the incident flux
    !! up is defined as negative in the adding method, so we make it positive here
    flx_down(:) = (D(:,1) + T(:))*mu_in(nlev)*F0_in
    flx_up(:) = -U(:,1)*mu_in(nlev)*F0_in

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