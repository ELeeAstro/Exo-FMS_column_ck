!!!
! Elspeth KH Lee - Aug 2023 : Initial version
! lw: Extended Absorption Approximation (EAA) following Li (2002) - This is the linear in tau method
!     Pros: Very fast method with LW scattering approximation, no matrix inversions
!     Cons: Not technically multiple scattering (can be quite innacurate at even moderate albedo)
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!
!!!

module lw_AA_L_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4
  use BB_integrate_mod, only : BB_integrate
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: sb = 5.670374419e-8_dp

  !! Gauss quadrature variables, cosine angle values (uarr) and weights (w)
  !! here you can comment in/out groups of mu values for testing
  !! make sure to make clean and recompile if you change these.

  !! Optimised quadrature for 1 node - 2 stream (Hogan 2024)
  ! integer, parameter :: nmu = 1
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.6096748751_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/1.0_dp/)

  !! Gauss–Jacobi-5 quadrature for 2 nodes - 4 stream (Hogan 2024)
  !integer, parameter :: nmu = 2
  !real(dp), dimension(nmu), parameter :: uarr = (/0.2509907356_dp, 0.7908473981_dp/)
  !real(dp), dimension(nmu), parameter :: w = (/0.2300253764_dp, 0.7699746236_dp/)

  !! Gauss–Jacobi-5 quadrature for 3 nodes - 6 stream (Hogan 2024)
  ! integer, parameter :: nmu = 3
  ! real(dp), dimension(nmu), parameter :: uarr = (/0.1024922169_dp, 0.4417960320_dp, 0.8633751621_dp/)
  ! real(dp), dimension(nmu), parameter :: w = (/0.0437820218_dp, 0.3875796738_dp, 0.5686383044_dp/)

  !! Gauss–Jacobi-5 quadrature for 4 nodes - 8 stream (Hogan 2024)
  integer, parameter :: nmu = 4
  real(dp), dimension(nmu), parameter :: uarr = (/0.0454586727_dp, 0.2322334416_dp, 0.5740198775_dp, 0.9030775973_dp/)
  real(dp), dimension(nmu), parameter :: w = (/0.0092068785_dp, 0.1285704278_dp, 0.4323381850_dp, 0.4298845087_dp /)

  private :: lw_AA_linear
  public :: lw_AA_L

contains


  subroutine lw_AA_L(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
        call lw_AA_linear(nlay, nlev, be(b,:), be_int(b), .false., tau_e(g,b,:), ssa(g,b,:), gg(g,b,:), nmu, uarr, w, &
          & lw_up_g(:), lw_down_g(:))
        lw_up(:) = lw_up(:) + lw_up_g(:) * gw(g)
        lw_down(:) = lw_down(:) + lw_down_g(:) * gw(g)
      end do
    end do

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_AA_L

  subroutine lw_AA_linear(nlay, nlev, be, be_int, surf, tau_in, w_in, g_in, nm, mu_q, w_q, &
    & flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev, nm
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int
    logical, intent(in) :: surf                     
    real(dp), dimension(nm), intent(in) :: mu_q, w_q 

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables
    integer :: k, m
    real(dp), dimension(nlay) :: dtau, eps, dtau_a
    real(dp), dimension(nlay) :: del, T, e0i, psi_u, psi_0
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g
    real(dp), parameter :: small_tau = 1.0e-3_dp

    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)

    !! EAA modified co-albedo.  No delta scaling is needed because this
    !! transport combination is invariant under a consistent delta transform.
    eps(:) = sqrt((1.0_dp - w_in(:))*(1.0_dp - g_in(:)*w_in(:)))

    !! Absorption/modified optical depth for transmission function
    dtau_a(:) = eps(:)*dtau(:)

    ! Zero the total flux arrays
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nm

      !! Use their Taylor expansions at small slant absorption optical thickness.
      del(:) = dtau_a(:)/mu_q(m)
      where (del(:) < small_tau)
        psi_0(:) = del(:)*(0.5_dp + del(:)*(-1.0_dp/6.0_dp + del(:)/24.0_dp))
        psi_u(:) = del(:)*(0.5_dp + del(:)*(-1.0_dp/3.0_dp + del(:)/8.0_dp))
        T(:) = 1.0_dp - psi_u(:) - psi_0(:)
      elsewhere
        T(:) = exp(-del(:))
        e0i(:) = 1.0_dp - T(:)
        psi_0(:) = (del(:) - e0i(:))/del(:)
        psi_u(:) = e0i(:) - psi_0(:)
      end where

      !! Begin angular sweeps
      !! Perform downward loop first
      ! Top boundary condition - no downward intensity at the top
      lw_down_g(1) = 0.0_dp
      do k = 1, nlay
        lw_down_g(k+1) = lw_down_g(k)*T(k) + psi_u(k)*be(k) + psi_0(k)*be(k+1)
      end do

      !! Perform upward loop
      ! Lower boundary condition - internal heat definition Fint = F_up - F_down,
      ! or a black-body surface (surf)
      if (surf) then
        lw_up_g(nlev) = be_int
      else
        lw_up_g(nlev) = lw_down_g(nlev) + be_int
      end if
      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*T(k) + psi_0(k)*be(k) + psi_u(k)*be(k+1)
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      flx_down(:) = flx_down(:) + lw_down_g(:) *  w_q(m)
      flx_up(:) = flx_up(:) + lw_up_g(:) *  w_q(m)

    end do

    !! The flux is the integrated intensity * pi (in this GJ weighting scheme)
    flx_down(:) = pi * flx_down(:)
    flx_up(:) = pi * flx_up(:)

  end subroutine lw_AA_linear

end module lw_AA_L_mod
