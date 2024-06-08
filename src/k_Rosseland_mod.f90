!!!
! E.K.H. Lee - Sep 2020
! A module that helps with grey and non-grey opacities from the literature
! 1. Tan & Komacek (2019)'s UHJ grey parameter fit
! 2. Freedman et al. (2014)'s & Valencia et al. (2013)'s fitting function
! 3. Parmentier and Guillot (2014,2015) 3V band, 2IR picket fence scheme parameters, and Bond albedo function
! Can be used in 1D RCE/RT or 3D GCM codes
!! TODO: Interpolate Freedman et al. (2014) tables directly as an option
!!!

module k_Rosseland_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: onedivpi = 1.0_dp/pi
  real(dp), parameter :: deg_to_rad = pi/180.0_dp

  !! Coefficent parameters for Freedman et al. (2014) table fit
  real(dp), parameter :: c1 = 10.602_dp
  real(dp), parameter :: c2 = 2.882_dp
  real(dp), parameter :: c3 = 6.09e-15_dp
  real(dp), parameter :: c4 = 2.954_dp
  real(dp), parameter :: c5 = -2.526_dp
  real(dp), parameter :: c6 = 0.843_dp
  real(dp), parameter :: c7 = -5.490_dp
  real(dp), parameter :: c8_l = -14.051_dp, c8_h = 82.241_dp
  real(dp), parameter :: c9_l = 3.055_dp, c9_h = -55.456_dp
  real(dp), parameter :: c10_l = 0.024_dp, c10_h = 8.754_dp
  real(dp), parameter :: c11_l = 1.877_dp, c11_h = 0.7048_dp
  real(dp), parameter :: c12_l = -0.445_dp, c12_h = -0.0414_dp
  real(dp), parameter :: c13_l = 0.8321_dp, c13_h = 0.8321_dp

  !! Coefficents parameters for the Valencia et al. (2013) table fit
  real(dp), parameter :: c1_v = -37.50_dp
  real(dp), parameter :: c2_v = 0.00105_dp
  real(dp), parameter :: c3_v = 3.2610_dp
  real(dp), parameter :: c4_v = 0.84315_dp
  real(dp), parameter :: c5_v = -2.339_dp
  real(dp), parameter :: c6_vl = -14.051_dp, c6_vh = 82.241_dp
  real(dp), parameter :: c7_vl = 3.055_dp, c7_vh = -55.456_dp
  real(dp), parameter :: c8_vl = 0.024_dp, c8_vh = 8.754_dp
  real(dp), parameter :: c9_vl = 1.877_dp, c9_vh = 0.7048_dp
  real(dp), parameter :: c10_vl = -0.445_dp, c10_vh = -0.0414_dp
  real(dp), parameter :: c11_vl = 0.8321_dp, c11_vh = 0.8321_dp

  !private
  public :: k_Ross_Freedman, gam_Parmentier, Bond_Parmentier, k_Ross_TK19, k_Ross_Valencia

contains

  !! Calculates the pressure dependent V and IR Rosseland mean opacity acording to the
  !! Tan & Komacek (2019) UHJ methods
  subroutine k_Ross_TK19(P, k_V, k_IR)
    implicit none

    ! Input:
    ! P - gas pressure [Pa]

    ! Output:
    ! k_V - V band Rosseland mean opacity [m2 kg-1]
    ! k_IR - IR band Rosseland mean pacity [m2 kg-1]

    real(dp), intent(in) :: P
    real(dp), intent(out) :: k_V, k_IR

    real(dp) :: Pl10

    Pl10 = log10(P)

    k_V = 10.0_dp**(0.0478_dp*Pl10**2 - 0.1366_dp*Pl10 - 3.2095_dp)
    k_IR = 10.0_dp**(0.0498_dp*Pl10**2 - 0.1329_dp*Pl10 - 2.9457_dp)


  end subroutine k_Ross_TK19


  !! Calculates the IR band Rosseland mean opacity (local T) according to the
  !! Freedman et al. (2014) fit and coefficents
  subroutine k_Ross_Freedman(Tin, Pin, met, k_IR)
    implicit none

    ! Input:
    ! T - Local gas temperature [K]
    ! P - Local gas pressure [pa]
    ! met - Local metallicity [M/H] (log10 from solar, solar [M/H] = 0.0)

    ! Output:
    ! k_IR - IR band Rosseland mean opacity [m2 kg-1]

    real(dp), intent(in) :: Tin, Pin, met
    real(dp), intent(out) :: k_IR

    real(dp) :: k_lowP, k_hiP
    real(dp) :: T, P, Tl10, Pl10

    T = Tin
    P = Pin * 10.0_dp ! Convert to dyne cm-2

    Tl10 = log10(T)
    Pl10 = log10(P)

    ! Low pressure expression
    k_lowP = c1*atan(Tl10 - c2) &
      & - (c3/(Pl10 + c4))*exp((Tl10 - c5)**2) &
      & + c6*met + c7

    ! Temperature split for coefficents = 800 K
    if (T <= 800.0_dp) then
      k_hiP = c8_l + c9_l*Tl10 &
        & + c10_l*Tl10**2 + Pl10*(c11_l + c12_l*Tl10) &
        & + c13_l * met * (0.5_dp + onedivpi*atan((Tl10 - 2.5_dp) / 0.2_dp))
    else
      k_hiP = c8_h + c9_h*Tl10 &
        & + c10_h*Tl10**2 + Pl10*(c11_h + c12_h*Tl10) &
        & + c13_h * met * (0.5_dp + onedivpi*atan((Tl10 - 2.5_dp) / 0.2_dp))
    end if

    ! Total Rosseland mean opacity - converted to m2 kg-1
    k_IR = (10.0_dp**k_lowP + 10.0_dp**k_hiP) / 10.0_dp

    ! Avoid divergence in fit for large values
    k_IR = min(k_IR,1.0e30_dp)

  end subroutine k_Ross_Freedman

  !! Calculates the IR band Rosseland mean opacity (local T) according to the
  !! Valencia et al. (2013) fit and coefficents
  subroutine k_Ross_Valencia(Tin, Pin, met, k_IR)
    implicit none

    ! Input:
    ! T - Local gas temperature [K]
    ! P - Local gas pressure [pa]
    ! met - Local metallicity [M/H] (log10 from solar, solar [M/H] = 0.0)

    ! Output:
    ! k_IR - IR band Rosseland mean opacity [m2 kg-1]

    real(dp), intent(in) :: Tin, Pin, met
    real(dp), intent(out) :: k_IR

    real(dp) :: k_lowP, k_hiP
    real(dp) :: Tl10, Pl10
    real(dp) :: T, P


    T = Tin
    P = Pin * 10.0_dp ! Convert to dyne

    Tl10 = log10(T)
    Pl10 = log10(P)

    k_lowP = c1_v * (Tl10-c2_v*Pl10-c3_v)**2 + (c4_v*met + c5_v)

    ! Temperature split for coefficents = 800 K
    if (T <= 800.0_dp) then
      k_hiP = (c6_vl+c7_vl*Tl10+c8_vl*Tl10**2) &
      &     + Pl10*(c9_vl+c10_vl*Tl10) &
      &     + met*c11_vl*(0.5_dp + onedivpi*atan((Tl10-2.5_dp)/0.2_dp))
    else
      k_hiP = (c6_vh+c7_vh*Tl10+c8_vh*Tl10**2) &
      &     + Pl10*(c9_vh+c10_vh*Tl10) &
      &     + met*c11_vh*(0.5_dp + onedivpi*atan((Tl10-2.5_dp)/0.2_dp))

    end if

    ! Total Rosseland mean opacity - converted to m2 kg-1
    k_IR = (10.0_dp**k_lowP + 10.0_dp**k_hiP) / 10.0_dp

    ! Avoid divergence in fit for large values
    k_IR = min(k_IR,1.0e30_dp)

  end subroutine k_Ross_Valencia

  !! Calculates 3 band grey visual gamma values and 2 picket fence IR gamma values
  !! according to the coefficents and equations in:
  !! Parmentier & Menou (2014) and Parmentier et al. (2015)
  !! NOTE: This does not calculate the opacity - call k_Ross_Freedman for that
  subroutine gam_Parmentier(Teff, table_num, gam_V, Beta_V, Beta, gam_1, gam_2, gam_P, tau_lim)
    implicit none

    ! Input:
    ! Teff - Effective temperature [K] (See Parmentier papers for various ways to calculate this)
    ! for non-irradiated atmosphere Teff = Tint
    ! table_num - Table selection from Parmentier et al. (2015): 1 = w. TiO/VO, 2 = w.o. TiO/VO

    ! Output:
    ! gam_V(3) - gamma ratio for 3 visual bands (gam_V = kV_Ross/kIR_Ross)
    ! beta_V(3) - fraction of total incident stellar flux in band (1/3 for Parmentier values)
    ! Beta - equilvalent bandwidth for picket fence IR model
    ! gam_1 - gamma ratio for IR band 1 (gam_1 = kIR_1/kIR_Ross)
    ! gam_2 - gamma ratio for IR band 2 (gam_2 = kIR_2/kIR_Ross)
    ! gam_P - gamma ratio for Planck mean (gam_P = kIR_Planck/kIR_Ross)
    ! tau_lim - tau limit variable (usually for IC system)

    real(dp), intent(in) :: Teff
    integer, intent(in) :: table_num

    real(dp), dimension(3), intent(out) :: gam_V, Beta_V
    real(dp), dimension(2), intent(out) :: Beta
    real(dp), intent(out) :: gam_1, gam_2, gam_P, tau_lim

    real(dp) :: R
    real(dp) :: aP, bP, cP
    real(dp) :: aV1, bV1, aV2, bV2, aV3, bV3
    real(dp) :: aB, bB
    real(dp) :: l10T, l10T2, RT

    ! Log 10 T_eff variables
    l10T = log10(Teff)
    l10T2 = l10T**2

    if (table_num == 1) then
      ! First table in Parmentier et al. (2015) w. TiO/VO
      ! Start large if statements with visual band and Beta coefficents
      if (Teff <= 200.0_dp) then
        aV1 = -5.51_dp ; bV1 = 2.48_dp
        aV2 = -7.37_dp ; bV2 = 2.53_dp
        aV3 = -3.03_dp ; bV3 = -0.20_dp
        aB = 0.84_dp  ; bB = 0.0_dp
      else if ((Teff > 200.0_dp) .and. (Teff <= 300.0_dp)) then
        aV1 = 1.23_dp ; bV1 = -0.45_dp
        aV2 = 13.99_dp ; bV2 = -6.75_dp
        aV3 = -13.87_dp ; bV3 = 4.51_dp
        aB = 0.84_dp  ; bB = 0.0_dp
      else if ((Teff > 300.0_dp) .and. (Teff <= 600.0_dp)) then
        aV1 = 8.65_dp ; bV1 = -3.45_dp
        aV2 = -15.18_dp ; bV2 = 5.02_dp
        aV3 = -11.95_dp ; bV3 = 3.74_dp
        aB = 0.84_dp  ; bB = 0.0_dp
      else if ((Teff > 600.0_dp) .and. (Teff <= 1400.0_dp)) then
        aV1 = -12.96_dp ; bV1 = 4.33_dp
        aV2 = -10.41_dp ; bV2 = 3.31_dp
        aV3 = -6.97_dp ; bV3 = 1.94_dp
        aB = 0.84_dp  ; bB = 0.0_dp
      else if ((Teff > 1400.0_dp) .and. (Teff < 2000.0_dp)) then
        aV1 = -23.75_dp ; bV1 = 7.76_dp
        aV2 = -19.95_dp ; bV2 = 6.34_dp
        aV3 = -3.65_dp ; bV3 = 0.89_dp
        aB = 0.84_dp  ; bB = 0.0_dp
      else if (Teff >= 2000.0_dp) then
        aV1 = 12.65_dp ; bV1 = -3.27_dp
        aV2 = 13.56_dp ; bV2 = -3.81_dp
        aV3 = -6.02_dp ; bV3 = 1.61_dp
        aB = 6.21_dp  ; bB = -1.63_dp
      end if

      ! gam_P coefficents
      aP = -2.36_dp
      bP = 13.92_dp
      cP = -19.38_dp

    else if (table_num == 2) then
      ! Appendix table from Parmentier et al. (2015) - without TiO and VO
      if (Teff <= 200.0_dp) then
        aV1 = -5.51_dp ; bV1 = 2.48_dp
        aV2 = -7.37_dp ; bV2 = 2.53_dp
        aV3 = -3.03_dp ; bV3 = -0.20_dp
        aB = 0.84_dp  ; bB = 0.0_dp
      else if ((Teff > 200.0_dp) .and. (Teff <= 300.0_dp)) then
        aV1 = 1.23_dp ; bV1 = -0.45_dp
        aV2 = 13.99_dp ; bV2 = -6.75_dp
        aV3 = -13.87_dp ; bV3 = 4.51_dp
        aB = 0.84_dp  ; bB = 0.0_dp
      else if ((Teff > 300.0_dp) .and. (Teff <= 600.0_dp)) then
        aV1 = 8.65_dp ; bV1 = -3.45_dp
        aV2 = -15.18_dp ; bV2 = 5.02_dp
        aV3 = -11.95_dp ; bV3 = 3.74_dp
        aB = 0.84_dp  ; bB = 0.0_dp
      else if ((Teff > 600.0_dp) .and. (Teff <= 1400.0_dp)) then
        aV1 = -12.96_dp ; bV1 = 4.33_dp
        aV2 = -10.41_dp ; bV2 = 3.31_dp
        aV3 = -6.97_dp ; bV3 = 1.94_dp
        aB = 0.84_dp  ; bB = 0.0_dp
      else if ((Teff > 1400.0_dp) .and. (Teff < 2000.0_dp)) then
        aV1 = -1.68_dp ; bV1 = 0.75_dp
        aV2 = 6.96_dp ; bV2 = -2.21_dp
        aV3 = 0.02_dp ; bV3 = -0.28_dp
        aB = 3.0_dp  ; bB = -0.69_dp
      else if (Teff >= 2000.0_dp) then
        aV1 = 10.37_dp ; bV1 = -2.91_dp
        aV2 = -2.4_dp ; bV2 = 0.62_dp
        aV3 = -16.54_dp ; bV3 = 4.74_dp
        aB = 3.0_dp  ; bB = -0.69_dp
      end if

      !gam_P coefficents
      if (Teff <= 1400.0_dp) then
        aP = -2.36_dp
        bP = 13.92_dp
        cP = -19.38_dp
      else
        aP = -12.45_dp
        bP = 82.25_dp
        cP = -134.42_dp
      end if

    end if

    ! Calculation of all values
    ! Visual band gamma
    gam_V(1) = 10.0_dp**(aV1 + bV1 * l10T)
    gam_V(2) = 10.0_dp**(aV2 + bV2 * l10T)
    gam_V(3) = 10.0_dp**(aV3 + bV3 * l10T)

    ! Visual band fractions
    Beta_V(:) = 1.0_dp/3.0_dp

    ! gamma_Planck - if < 1 then make it grey approximation (k_Planck = k_Ross, gam_P = 1)
    gam_P = 10.0_dp**(aP * l10T2 + bP * l10T + cP)
    if (gam_P < 1.0000001_dp) then
      gam_P = 1.0000001_dp
    end if

    ! equivalent bandwidth value
    Beta(1) = aB + bB * l10T
    Beta(2) = 1.0_dp - Beta(1)

    ! IR band kappa1/kappa2 ratio - Eq. 96 from Parmentier & Menou (2014)
    RT = (gam_P - 1.0_dp)/(2.0_dp*Beta(1)*Beta(2))
    R = 1.0_dp + RT + sqrt(RT**2 + RT)

    ! gam_1 and gam_2 values - Eq. 92, 93 from Parmentier & Menou (2014)
    gam_1 = Beta(1) + R - Beta(1)*R
    gam_2 = gam_1 / R

    ! Calculate tau_lim parameter
    tau_lim = 1.0_dp/(gam_1*gam_2) * sqrt(gam_P/3.0_dp)

  end subroutine gam_Parmentier

  ! Calculates the Bond Albedo according to Parmentier et al. (2015) expression
  subroutine Bond_Parmentier(Teff0, grav,  AB)
    implicit none

    ! Input:
    ! Teff0 - Atmospheric profile effective temperature [K] with zero albedo
    ! grav - Surface gravity of planet [m s-2]

    ! Output:
    ! AB - Bond albedo

    real(dp), intent(in) :: Teff0, grav
    real(dp), intent(out) :: AB

    real(dp) :: a, b

    ! a and b cofficents dependent on T_eff and grav
    if (Teff0 <= 250.0_dp) then
      a = -0.335_dp * grav**(0.070_dp)
      b = 0.0_dp
    else if ((Teff0 > 250.0_dp) .and. (Teff0 <= 750.0_dp)) then
      a = -0.335_dp * grav**(0.070_dp) + 2.149_dp * grav**(0.135_dp)
      b = -0.896_dp * grav**(0.135_dp)
    else if ((Teff0 > 750.0_dp) .and. (Teff0 < 1250.0_dp)) then
      a = -0.335_dp * grav**(0.070_dp) -  0.428_dp * grav**(0.135_dp)
      b = 0.0_dp
    else if (Teff0 >= 1250.0_dp) then
      a = 16.947_dp - 3.174_dp * grav**(0.070_dp) - 4.051_dp * grav**(0.135_dp)
      b = -5.472_dp + 0.917_dp * grav**(0.070_dp) + 1.170_dp * grav**(0.135_dp)
    end if

    ! Final Bond Albedo expression
    AB = 10.0_dp**(a + b * log10(Teff0))

  end subroutine Bond_Parmentier

end module k_Rosseland_mod

!! A test program for the above module to show example call stuctures
!! Comment out when coupling to other codes
! program test_k_Rosseland_mod
!   use k_Rosseland_mod, only: k_Ross_Freedman, gam_Parmentier, Bond_Parmentier, k_Ross_Tan, k_Ross_Valencia
!   implicit none
!
!   double precision :: T = 1000.0d0
!   double precision :: P = 1.0d0
!   double precision :: met = 0.0d0
!
!   double precision :: k_V_Tan, k_IR_Tan, k_IR_Freedman, k_IR_Valencia
!
!   call k_Ross_Tan(P*1d5, k_V_Tan, k_IR_Tan)
!   print*, 'Tan: ', T, P, k_V_Tan, k_IR_Tan
!
!   call k_Ross_Freedman(T, P*1d5, met, k_IR_Freedman)
!   print*, 'Freedman: ', T, P, met, k_IR_Freedman
!
!   call k_Ross_Valencia(T, P*1d5, met, k_IR_Valencia)
!   print*, 'Valencia: ', T, P, met, k_IR_Valencia
!
! end program test_k_Rosseland_mod
