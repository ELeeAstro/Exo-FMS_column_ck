!!!
! Elspeth KH Lee - Aug 2023 : Initial version
!                - Nov 2026 : Converted to corr-k (band + g-ordinance loops) from Exo-FMS_column_semigrey
!                - Nov 2026 : Optimisations - log(B ratio) once per layer (not per angle and sweep),
!                             reciprocal angles, contiguous per g-ordinance inputs
!
! lw: (Extended) Absorption Approximation following Li (2002) (alpha-nEAA) - This is the exoponential in tau method
!     Pros: Very fast method with LW scattering approximation, no matrix inversions
!     Cons: Not technically multiple scattering (can be quite innacurate at even moderate albedo)
!
!!!

module lw_AA_E_mod
  use, intrinsic :: iso_fortran_env
  use WENO4_mod, only : interpolate_weno4
  use BB_integrate_mod, only : BB_integrate
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Required constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: sb = 5.670374419e-8_dp
  real(dp), parameter :: dtau_floor = 1.0e-12_dp

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
  real(dp), dimension(nmu), parameter :: ruarr = 1.0_dp/uarr

  private :: lw_AA_exp, aa_source_down, aa_source_up, exprel
  public :: lw_AA_E

contains

  subroutine lw_AA_E(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
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
        ! Contiguous copies (avoids array temporaries for the strided slices)
        tau_g(:) = tau_e(g,b,:)
        ssa_g(:) = ssa(g,b,:)
        gg_g(:) = gg(g,b,:)
        call lw_AA_exp(nlay, nlev, be_b, be_int(b), tau_g, ssa_g, gg_g, lw_up_g(:), lw_down_g(:))
        lw_up(:) = lw_up(:) + lw_up_g(:) * gw(g)
        lw_down(:) = lw_down(:) + lw_down_g(:) * gw(g)
      end do
    end do

    !! Net lw fluxes at each level
    lw_net(:) = lw_up(:) - lw_down(:)

    !! Outgoing Longwave Radiation (olr)
    olr = lw_up(1)

  end subroutine lw_AA_E

  subroutine lw_AA_exp(nlay, nlev, be, be_int, tau_in, w_in, g_in, flx_up, flx_down)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlev), intent(in) :: be, tau_in
    real(dp), dimension(nlay), intent(in) :: w_in, g_in
    real(dp), intent(in) :: be_int

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: flx_up, flx_down

    !! Work variables
    integer :: k, m
    real(dp), dimension(nlay) :: dtau, eps, dtau_a, lr
    real(dp), dimension(nlay) :: T, x
    real(dp), dimension(nlev) :: lw_up_g, lw_down_g

    dtau(:) = tau_in(2:nlev) - tau_in(1:nlay)

    !! EAA modified co-albedo.  No delta scaling is needed because this
    !! transport combination is invariant under a consistent delta transform.
    eps(:) = sqrt((1.0_dp - w_in(:))*(1.0_dp - g_in(:)*w_in(:)))

    !! Absorption/modified optical depth for transmission function
    dtau_a(:) = eps(:)*dtau(:)

    !! log of the Planck ratio across each layer (angle independent)
    lr(:) = log(be(2:nlev)/be(1:nlay))

    ! Zero the total flux arrays
    flx_up(:) = 0.0_dp
    flx_down(:) = 0.0_dp

    !! Start loops to integrate in mu space
    do m = 1, nmu

      !! Efficency variables
      x(:) = dtau_a(:)*ruarr(m)
      T(:) = exp(-x(:))

      !! Begin two-stream loops
      !! Perform downward loop first - also calculate efficency variables
      ! Top boundary condition - 0 flux downward from top boundary
      lw_down_g(1) = 0.0_dp
      do k = 1, nlay
        lw_down_g(k+1) = lw_down_g(k)*T(k) + &
          & aa_source_down(be(k), be(k+1), x(k), T(k), lr(k))
      end do

      !! Perform upward loop
      ! Gas-giant lower boundary: F_up - F_down = sigma*Tint**4.
      lw_up_g(nlev) = lw_down_g(nlev) + be_int

      do k = nlay, 1, -1
        lw_up_g(k) = lw_up_g(k+1)*T(k) + &
          & aa_source_up(be(k), be(k+1), x(k), T(k), lr(k))
      end do

      !! Sum up flux arrays with Gaussian quadrature weights and points for this mu stream
      flx_down(:) = flx_down(:) + lw_down_g(:) * w(m)
      flx_up(:) = flx_up(:) + lw_up_g(:) * w(m)

    end do

    !! The flux is the integrated intensity * pi (in this GJ weighting scheme)
    flx_down(:) = pi * flx_down(:)
    flx_up(:) = pi * flx_up(:)
  

  end subroutine lw_AA_exp

  !! log_ratio = log(bbot/btop)
  pure elemental function aa_source_down(btop, bbot, x, trans, log_ratio) result(source)
    implicit none

    real(dp), intent(in) :: btop, bbot, x, trans, log_ratio
    real(dp) :: source
    real(dp) :: z

    z = x + log_ratio

    if (abs(z) < 1.0e-4_dp) then
      source = btop*trans*x*exprel(z)
    else
      source = (x/z)*(bbot - btop*trans)
    end if

  end function aa_source_down

  !! log_ratio = log(bbot/btop)
  pure elemental function aa_source_up(btop, bbot, x, trans, log_ratio) result(source)
    implicit none

    real(dp), intent(in) :: btop, bbot, x, trans, log_ratio
    real(dp) :: source
    real(dp) :: z

    z = log_ratio - x

    if (abs(z) < 1.0e-4_dp) then
      source = btop*x*exprel(z)
    else
      source = (x/(x - log_ratio))*(btop - bbot*trans)
    end if

  end function aa_source_up

  pure elemental function exprel(x) result(value)
    implicit none

    real(dp), intent(in) :: x
    real(dp) :: value

    if (abs(x) < 1.0e-4_dp) then
      value = 1.0_dp + x*(0.5_dp + x*(1.0_dp/6.0_dp + &
        & x*(1.0_dp/24.0_dp + x/120.0_dp)))
    else
      value = (exp(x) - 1.0_dp)/x
    end if

  end function exprel

end module lw_AA_E_mod
