!!!
! Elspeth KH Lee - May 2021
! Two-stream DISORT version, modified by Xianyu Tan to include an internal heat source.
! Pros: Stable, performs accurate scattering calculations tried and tested, reliable model.
! Cons: Slower than other methods.
!!!

module ts_disort_scatter_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  public :: ts_disort_scatter
  private :: linear_log_interp

contains

  subroutine ts_disort_scatter(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, mu_z, Finc, Tint, &
    & olr, net_F)
    implicit none

    !! Input variables
    integer, intent(in) :: nlay, nlev, nb, ng
    real(dp), dimension(ng), intent(in) :: gw
    real(dp), dimension(nb+1), intent(in) :: wn_e
    real(dp), dimension(nlay), intent(in) :: Tl, pl
    real(dp), dimension(nlev), intent(in) :: pe
    real(dp), dimension(ng,nb,nlev), intent(in) :: tau_e
    real(dp), dimension(ng,nb,nlay), intent(in) :: ssa, gg
    real(dp), dimension(nb), intent(in) :: Finc
    real(dp), intent(in) :: mu_z, Tint

    !! Output variables
    real(dp), dimension(nlev), intent(out) :: net_F
    real(dp), intent(out) :: olr

    !! Work variables
    integer :: i, b, g
    real(dp), dimension(nlev) :: Te
    real(dp) :: a_surf = 0.0_dp
    logical :: Bezier = .True.

    !! Conversion arrays from FMS to DISORT dimensions and work variables
    integer, parameter :: maxcly=200, maxulv=201
    real(dp), dimension(0:maxcly) :: Te_0
    real(dp) :: wvnmlo, wvnmhi
    real(dp), dimension(maxcly) :: dtauc, utau
    real(dp), dimension(maxcly) :: ggg, ssalb
    real(dp), dimension(maxulv) :: sw_net, lw_net
    real(dp), dimension(nb,maxulv) :: sw_net_b, lw_net_b
    real(dp), dimension(ng,maxulv) :: sw_net_g, lw_net_g
    real(dp) :: umu0, fbeam
    logical :: planck

    !! Find temperature at layer edges through interpolation and extrapolation
    if (Bezier .eqv. .True.) then
      ! Perform interpolation using Bezier peicewise polynomial interpolation
      do i = 2, nlay-1
        call bezier_interp(pl(i-1:i+1), Tl(i-1:i+1), 3, pe(i), Te(i))
        !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
      end do
      call bezier_interp(pl(nlay-2:nlay), Tl(nlay-2:nlay), 3, pe(nlay), Te(nlay))
    else
      ! Perform interpolation using linear interpolation
      do i = 2, nlay
        call linear_log_interp(pe(i), pl(i-1), pl(i), Tl(i-1), Tl(i), Te(i))
        !print*, i, pl(i), pl(i-1), pe(i), Tl(i-1), Tl(i), Te(i)
      end do
    end if

    ! Edges are linearly interpolated
    Te(1) = 10.0_dp**(log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))
    Te(nlev) = 10.0_dp**(log10(Tl(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))

    Te_0(0:nlay) = Te(1:nlev)

    !! Shortwave flux calculation
    if (mu_z > 0.0_dp) then
      planck = .False.
      umu0 = mu_z
      sw_net(:) = 0.0_dp
      do b = 1, nb
        sw_net_b(b,:) = 0.0_dp
        fbeam = Finc(b)
        wvnmlo = wn_e(b+1)
        wvnmhi = wn_e(b)
        do g = 1, ng
          ggg(1:nlay) = gg(g,b,:)
          ssalb(1:nlay) = ssa(g,b,:)
          utau(1:nlev) = tau_e(g,b,:)
          do i = 1, nlay
            dtauc(i) = (tau_e(g,b,i+1) - tau_e(g,b,i))
          end do
          call CALL_TWOSTR (nlay,Te_0,ggg,ssalb,dtauc,nlev,utau,planck,wvnmlo,wvnmhi,Tint,fbeam,umu0,sw_net_g(g,:))
          sw_net_b(b,:) = sw_net_b(b,:) + sw_net_g(g,:) * gw(g)
        end do
        sw_net(:) = sw_net(:) + sw_net_b(b,:)
      end do
    else
      sw_net(:) = 0.0_dp
    end if

    !! Longwave two-stream flux calculation
    planck = .True.
    fbeam = 0.0_dp
    umu0 = 1.0_dp
    lw_net(:) = 0.0_dp
    do b = 1, nb
      lw_net_b(b,:) = 0.0_dp
      wvnmlo = wn_e(b+1)
      wvnmhi = wn_e(b)
      do g = 1, ng
        ggg(1:nlay) = gg(g,b,:)
        ssalb(1:nlay) = ssa(g,b,:)
        utau(1:nlev) = tau_e(g,b,:)
        do i = 1, nlay
          dtauc(i) = (tau_e(g,b,i+1) - tau_e(g,b,i))
        end do
        call CALL_TWOSTR (nlay,Te,ggg,ssalb,dtauc,nlev,utau,planck,wvnmlo,wvnmhi,Tint,fbeam,umu0,lw_net_g(g,:))
        lw_net_b(b,:) = lw_net_b(b,:) + lw_net_g(g,:) * gw(g)
      end do
      lw_net(:) = lw_net(:) + lw_net_b(b,:)
    end do

    !! Net fluxes at each level
    net_F(:) = lw_net(1:nlev) + sw_net(1:nlev)

    olr = 0.0_dp

  end subroutine ts_disort_scatter

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / log10(x2/x1)

    yval = 10.0_dp**((ly1 * log10(x2/xval) + ly2 * log10(xval/x1)) * norm)

  end subroutine linear_log_interp

  subroutine bezier_interp(xi, yi, ni, x, y)
    implicit none

    integer, intent(in) :: ni
    real(dp), dimension(ni), intent(in) :: xi, yi
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y

    real(dp) :: xc, dx, dx1, dy, dy1, w, yc, t, wlim, wlim1

    !xc = (xi(1) + xi(2))/2.0_dp ! Control point (no needed here, implicitly included)
    dx = xi(2) - xi(1)
    dx1 = xi(3) - xi(2)
    dy = yi(2) - yi(1)
    dy1 = yi(3) - yi(2)

    if (x > xi(1) .and. x < xi(2)) then
      ! left hand side interpolation
      !print*,'left'
      w = dx1/(dx + dx1)
      ! wlim = 1.0_dp + 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      ! wlim1 = 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      ! if (w < wlim .or. w > wlim1) then
      !  w = 1.0_dp
      ! end if
      yc = yi(2) - dx/2.0_dp * (w*dy/dx + (1.0_dp - w)*dy1/dx1)
      t = (x - xi(1))/dx
      y = (1.0_dp - t)**2 * yi(1) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(2)
    else ! (x > xi(2) and x < xi(3)) then
      ! right hand side interpolation
      !print*,'right'
      w = dx/(dx + dx1)
      ! wlim = 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      ! wlim1 = 1.0_dp + 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      ! if (w < wlim .or. w > wlim1) then
      !  w = 1.0_dp
      ! end if
      yc = yi(2) + dx1/2.0_dp * (w*dy1/dx1 + (1.0_dp - w)*dy/dx)
      t = (x - xi(2))/(dx1)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine bezier_interp

end module ts_disort_scatter_mod
