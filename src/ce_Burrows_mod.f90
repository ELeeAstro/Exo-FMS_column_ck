!!!
! Elspeth K.H. Lee (Jun 2024):
! Burrows et al. (1999) analytical chemical equilibrium method for 
! H2, He, H2O, CH4, CO, N2, NH3 
! Solar metallicity only!
!!!

module ce_Burrows_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  !! Coefficents for Burrows analytic expression
  real(dp), parameter :: a1 = 1.106131e6_dp, b1 = -5.6895e4_dp, c1 = 62.565_dp
  real(dp), parameter :: d1 = -5.81396e-4_dp, e1 = 2.346515e-8_dp
  real(dp), parameter :: a2 = 8.16413e5_dp, b2 = -2.9109e4_dp, c2 = 58.5878_dp
  real(dp), parameter :: d2 = -7.8284e-4_dp, e2 = 4.729048e-8_dp
  real(dp), parameter :: R_cal = 1.98720425864083_dp

contains

  subroutine analytic_Burrows(nsp, T_in, P_in, vmr_out)
    implicit none

    integer, intent(in) :: nsp
    real(dp), intent(in) :: T_in, P_in
    real(dp), dimension(nsp), intent(out) :: vmr_out

    real(dp) :: K1, K2, A_C, A_O, A_N, A_Si, P_H2
    real(dp) :: B_CO, B_CH4, B_H2O, B_N2, B_NH3
    real(dp) :: x_H2, x_He,  x_CO, x_CH4, x_H2O, x_N2, x_NH3


    ! Calculates 7 species in EQ analytically - solar only (for now)
    ! For this to work, species have to be in order in the array

    if (nsp /= 7) then
      print*, 'Error in CE_Burrows, nsp /= 7', nsp
      stop
    end if

    ! CO, CH4, H2O system
    K1 = exp((a1/T_in + b1 + c1*T_in + d1*T_in**2 + e1*T_in**3)/(R_cal * T_in))

    A_C = 10.0_dp**(8.50_dp - 12.0_dp)
    if (T_in < 1500.0_dp) then
      ! Remove a portion of O due to Si species condensation
      A_Si = 10.0_dp**(7.51_dp - 12.0_dp)
      A_O = 10.0_dp**(8.76_dp - 12.0_dp) - 3.28_dp * A_Si
    else
      ! No Si condensation, normal abundances
      A_O = 10.0_dp**(8.76_dp - 12.0_dp)
    end if

    P_H2 = 0.91183_dp * P_in
    x_H2 = 0.91183_dp

    B_CO = A_C + A_O + P_H2**2/(2.0_dp * K1) - sqrt((A_C + A_O + P_H2**2/(2.0_dp * K1))**2 - 4.0_dp*A_C*A_O)
    B_CH4 = 2.0_dp * A_C - B_CO
    B_H2O = 2.0_dp * A_O - B_CO

    x_CO = (B_CO * P_H2)/P_in; x_CH4 = (B_CH4 * P_H2)/P_in; x_H2O = (B_H2O * P_H2)/P_in

    ! N2, NH3 system
    K2 = exp((a2/T_in + b2 + c2*T_in + d2*T_in**2 + e2*T_in**3)/(R_cal * T_in))

    A_N = 10.0_dp**(7.86_dp - 12.0_dp)

    B_N2 = A_N + P_H2**2/(8.0_dp * K2) - sqrt((A_N + P_H2**2/(8.0_dp * K2))**2 - A_N**2)
    B_NH3 = 2.0_dp * (A_N - B_N2)

    x_N2 = (B_N2 * P_H2)/P_in
    x_NH3 = (B_NH3 * P_H2)/P_in

    ! Can assume Helium VMR is 1 - the other species VMR
    x_He = 1.0_dp - (x_H2 + x_H2O + x_CH4 + x_CO + x_N2 + x_NH3)

    vmr_out = (/x_H2, x_He, x_H2O, x_CH4, x_CO, x_N2, x_NH3/) 

  end subroutine analytic_Burrows

end module ce_Burrows_mod