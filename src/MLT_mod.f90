module MLT_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use WENO4_mod, only : interpolate_weno4  
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: alp = 1.0_dp, beta = 2.2_dp
  real(dp), parameter :: dt_max = 0.5_dp
  real(dp), parameter :: Kzz_min = 1e1_dp ! Minimum Kzz (typically 1e5 cm2 s-1)
  real(dp), parameter :: Kzz_max = 1e8_dp ! Maximum Kzz

  real(dp), parameter :: f1 = 1.0_dp/8.0_dp, f2 = 1.0_dp/2.0_dp

  public :: MLT

contains

  subroutine MLT(nlay, nlev, t_end, Tl, pl, pe, Rd_air, cp_air, kappa_air, &
    & grav, dT_conv, Kzz)
    implicit none

    integer, intent(in) :: nlay, nlev
    real(dp), dimension(nlay), intent(in) :: Rd_air, cp_air, kappa_air
    real(dp) ::  grav, t_end
    real(dp), dimension(nlay), intent(in) :: pl
    real(dp), dimension(nlev), intent(in) :: pe

    real(dp), dimension(nlay), intent(in) :: Tl
    real(dp), dimension(nlay), intent(out) :: dT_conv, Kzz

    integer :: k, ncall, krcb
    real(dp) :: Hp, L, rho, t_now, dt, w, w_over, w_rcb
    real(dp), dimension(nlay) :: gradrad, Fc, Tl_c
    real(dp), dimension(nlev) :: Fc_e, Te
    real(dp), dimension(nlay) :: Kzz_av, Kzz_ov

    !! Save the input temperature into copy
    Tl_c(:) = Tl(:)

    !! Loop over time until we reach t_end
    t_now = 0.0_dp
    ncall = 0
    Kzz_av(:) = 0.0_dp
    do while (t_now < t_end)

      !! Calculate the timestep
      if ((t_now + dt_max) >= t_end) then
        dt = t_end - t_now
      else
        dt = dt_max
      end if

      !! Use WENO4 method to (smoothly) interpolate layers to levels
      Te(:) = interpolate_weno4(pe, pl, Tl, .False.)

      !! Edges are linearly interpolated to avoid overshoot
      Te(1) = 10.0_dp**(log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))
      Te(nlev) = 10.0_dp**(log10(Tl(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * log10(Tl(nlay)/Te(nlay)))


      !! Calculate the regions that are convectivly unstable
      do k = nlay, 1, -1 
        gradrad(k) = log(Te(k)/Te(k+1))/log(pe(k)/pe(k+1))
      end do

      !! Follow the Joyce & Tayar (2023), Marley & Robinson (2015) and
      !Robinson & Marley (2014) formalisims
      do k = 1, nlay

        ! Scale height and char. length scale of layer
        Hp = (Rd_air(k) * Tl_c(k)) / grav
        L = alp * Hp

        if (gradrad(k) > kappa_air(k)) then
          ! Convective region, calculate convective flux

          ! Buoyancy characteristic vertical velocity Robinson & Marley (2014)
          w = (f1 * L**2 * grav * (gradrad(k) - kappa_air(k))) / Hp
          w = sqrt(w)
          ! Vertical convective thermal flux (Joyce & Tayar 2023)
          rho = pl(k)/(Rd_air(k) * Tl_c(k))
          Fc(k) = f2 * (rho * cp_air(k) * w * Tl_c(k) * L * (gradrad(k) - kappa_air(k)))/Hp

        else
          ! Convective flux is not adjusted
          Fc(k) = 0.0_dp
          w = 0.0_dp
        end if

        ! Keep running total of the Kzz value for this layer
        Kzz_av(k) = Kzz_av(k) + w * L        

      end do

      !! Use WENO4 method to (smoothly) interpolate layers to levels
      Fc_e(:) = interpolate_weno4(pe, pl, Fc, .False.)

      !! Edges are linearly extrapolted (or equal 0)
      Fc_e(1) = 0.0_dp !(Fc(1) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * (Fc(1) - Fc_e(2)))
      Fc_e(nlev) = 0.0_dp !(Fc(nlay)) + (log10(pe(nlev)/pe(nlay))/log10(pl(nlay)/pe(nlay))) * (Fc(nlay) - Fc_e(nlay))

      !! Calculate the temperature tendency due to radiation
      do k = 1, nlay
        dT_conv(k) = (grav/cp_air(k)) * (Fc_e(k+1)-Fc_e(k))/(pe(k+1) - pe(k))
        !print*, n, i, dT_conv(k)
      end do

      !! Forward march the temperature change copy   
      Tl_c(:) = Tl_c(:) + dt * dT_conv(:)

      !! Add dt to the time
      t_now = t_now + dt
      ncall = ncall + 1

    end do

    Kzz_av(:) = Kzz_av(:)/real(ncall,dp)

    !! Calculate the overshoot component - assume from end state (no running average)

    !! First find the rcb boundary 
    krcb = nlay
    w_rcb = 1e-30_dp
    do k = nlay, 1, -1
      if (Fc(k) > 0.0_dp) then
        krcb = k
      else
        krcb = k+1
        Hp = (Rd_air(krcb) * Tl_c(krcb)) / grav
        L = alp * Hp
        w_rcb = (f1 * L**2 * grav * (gradrad(krcb) - kappa_air(krcb))) / Hp
        w_rcb = sqrt(w_rcb)
        exit
      end if
    end do

    !! Now calculate the overshoot component
    Kzz_ov(:) = 0.0_dp
    do k = nlay, 1, -1
      if (k >= krcb) then
        !! In covective region - do not add overshoot
      else
        !! In overshoot region, add overshoot component
        w_over = exp(log(w_rcb) - beta * max(0.0_dp, log(pl(krcb)/pl(k))))
        Hp = (Rd_air(k) * Tl_c(k)) / grav
        L = alp * Hp
        Kzz_ov(k) = w_over * L
        if (Kzz_ov(k) < Kzz_min) then
          exit
        end if
      end if
    end do

    !! Make sure Kzz is above minimum value
    Kzz(:) = max(Kzz_av(:) + Kzz_ov(:),Kzz_min)

    !! Make sure Kzz is smaller than the maximum value
    Kzz(:) = min(Kzz(:),Kzz_max)

    !! Pass back the net temperature tendency
    dT_conv(:) = (Tl_c(:) - Tl(:))/t_end

  end subroutine MLT

end module MLT_mod
