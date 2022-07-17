!!!
! E.K.H. Lee - Sep 2020
!
!! TODO: Add boundary checks for pressure range
!!!

module CE_mod
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

  ! T-p and x data arrays from table
  integer :: nT, nP, npoint, nmol, u
  character(len=10), allocatable, dimension(:) :: m_name
  real(dp), allocatable, dimension(:) :: P_t, T_t, lP_t, lT_t
  real(dp), allocatable, dimension(:,:) :: mu_t
  real(dp), allocatable, dimension(:,:,:) :: x_t

  ! First call to module flag
  logical :: first_call = .True.

  public :: CE_interpolate_bilinear, CE_Burrows, CE_interpolate_Bezier
  private :: CE_interpolate_init, locate, linear_interp, bilinear_interp, Bezier_interp

contains

  subroutine CE_Burrows(T_in, P_in, m_size, x_out)
    implicit none

    integer, intent(in) :: m_size
    real(dp), intent(in) :: T_in, P_in
    real(dp), dimension(m_size), intent(out) :: x_out

    real(dp) :: K1, K2, A_C, A_O, A_N, A_Si, A_Ti, A_V, P_H2
    real(dp) :: B_CO, B_CH4, B_H2O, B_N2, B_NH3
    real(dp) :: x_H2, x_He,  x_CO, x_CH4, x_H2O, x_N2, x_NH3, x_TiO, x_VO


    ! Calculates 8 species in EQ analytically - solar only (for now)
    ! For this to work, species have to be in order in the array

    if (m_size /= 7) then
      print*, 'Error in CE_Burrows, m_size /= 8', m_size
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

!    if ((T_in > 1500.0_dp) .or. (P_in > 1e7_dp)) then
!      ! Simple Ti and V abundance and condensation scheme following Amundsen et al (2014)
!      A_Ti = 10.0_dp**(4.95_dp - 12.0_dp)
!      A_V = 10.0_dp**(3.93_dp - 12.0_dp)
!      x_TiO = A_Ti
!      x_VO = A_V
!    else
      x_TiO = 0.0_dp
      x_VO = 0.0_dp
!    end if

    ! Can assume Helium VMR is 1 - the other species VMR
    x_He = 1.0_dp - (x_H2 + x_H2O + x_CH4 + x_CO + x_N2 + x_NH3)

    x_out = (/x_H2, x_He, x_H2O,x_CH4,x_CO,x_N2,x_NH3/) !,x_TiO,x_VO/)

  end subroutine CE_Burrows

  subroutine CE_interpolate_bilinear(P_in, T_in, m_in, m_size, VMR_table, x_out, mu_out)
    implicit none

    !! Input:
    ! - m_size = size of molecule array
    ! - P_in = input pressure [bar]
    ! - T_in = input temperature [K]
    ! - m_in = name array of species to interpolate
    !! Output:
    ! - x_out = VMR of each species @ P_in and T_in
    ! - mu_out = mean molecular weight @ P_in and T_in
    integer, intent(in) :: m_size
    character(len=50) :: VMR_table
    real(dp), intent(in) :: P_in, T_in
    character(len=4), dimension(m_size), intent(in) :: m_in
    real(dp), intent(out), dimension(m_size) :: x_out
    real(dp), intent(out) :: mu_out

    !! Work variables
    integer :: m
    integer :: i_pl, i_pu, i_tl, i_tu
    real(dp) :: xval,x1,x2,yval,y1,y2,a11,a21,a12,a22,aval

    ! Initialise - read VMR from table for T and p grid
    if (first_call .eqv. .True.) then
      call CE_interpolate_init(VMR_table)
      first_call = .False.
    end if

   ! Find upper and lower T and P indexes
   call locate(P_t, nP, P_in, i_pl)
   i_pu = i_pl + 1
   call locate(T_t, nT, T_in, i_tl)
   i_tu = i_tl + 1

   ! print*, P_in
   ! print*, i_pl, i_pu, P_t(i_pl), P_t(i_pu)
   ! print*, T_t(i_tl), T_t(i_tu)
   ! print*, mu_t(i_pl,i_tu), mu_t(i_pu,i_tu)

   if (i_tl == nT) then
     !! Input higher than T grid boundary, make it = maxval(T)
     !! Perform mu linear interpolation
     xval = log10(P_in)
     x1 = lP_t(i_pl)
     x2 = lP_t(i_pu)
     y1 = mu_t(i_pl,nT)
     y2 = mu_t(i_pu,nT)
     call linear_interp(xval, x1, x2, y1, y2, yval)
     mu_out = yval

     do m = 1, m_size
       y1 = x_t(m,i_pl,nT)
       y2 = x_t(m,i_pu,nT)
       call linear_interp(xval, x1, x2, y1, y2, yval)
       x_out(m) = 10.0_dp**yval ! Unlog for output
     end do

   else if (i_tl == 0) then
     !! Input lower than T grid boundary, make it = minval(T)
     !! Perform mu linear interpolation
     xval = log10(P_in)
     x1 = lP_t(i_pl)
     x2 = lP_t(i_pu)
     y1 = mu_t(i_pl,1)
     y2 = mu_t(i_pu,1)
     call linear_interp(xval, x1, x2, y1, y2, yval)
     mu_out = yval

     do m = 1, m_size
       y1 = x_t(m,i_pl,1)
       y2 = x_t(m,i_pu,1)
       call linear_interp(xval, x1, x2, y1, y2, yval)
       x_out(m) = 10.0_dp**yval ! Unlog for output
     end do

   else
     !! Within T grid bounds
     !! Coordinates for T-p bi-linear interpolation
     xval = log10(T_in)
     x1 = lT_t(i_tl)
     x2 = lT_t(i_tu)
     yval = log10(P_in)
     y1 = lP_t(i_pl)
     y2 = lP_t(i_pu)

     !! Perform bi-linear interpolation for mu
     a11 = mu_t(i_pl,i_tl)
     a21 = mu_t(i_pu,i_tl)
     a12 = mu_t(i_pl,i_tu)
     a22 = mu_t(i_pu,i_tu)
     call bilinear_interp(xval, yval, x1, x2, y1, y2, a11, a21, a12, a22, aval)
     mu_out = aval

     !! Perform bi-linear interpolation for each species in a loop
     do m = 1, m_size
       a11 = x_t(m,i_pl,i_tl)
       a21 = x_t(m,i_pu,i_tl)
       a12 = x_t(m,i_pl,i_tu)
       a22 = x_t(m,i_pu,i_tu)
       call bilinear_interp(xval, yval, x1, x2, y1, y2, a11, a21, a12, a22, aval)
       x_out(m) = 10.0_dp**aval ! Unlog for output
     end do
   end if

  end subroutine CE_interpolate_bilinear


  subroutine CE_interpolate_Bezier(P_in, T_in, m_in, m_size, VMR_table, x_out, mu_out)
    implicit none

    !! Input:
    ! - m_size = size of molecule array
    ! - P_in = input pressure [bar]
    ! - T_in = input temperature [K]
    ! - m_in = name array of species to interpolate
    !! Output:
    ! - x_out = VMR of each species @ P_in and T_in
    ! - mu_out = mean molecular weight @ P_in and T_in

    integer, intent(in) :: m_size
    character(len=50) :: VMR_table
    real(dp), intent(in) :: P_in, T_in
    character(len=4), dimension(m_size), intent(in) :: m_in
    real(dp), intent(out), dimension(m_size) :: x_out
    real(dp), intent(out) :: mu_out

    integer :: i_p1, i_p2, i_p3, i_t1, i_t2, i_t3
    real(dp) :: lP_in, lT_in
    real(dp), dimension(3) :: lTa, lPa, mua, xa, xa_out, mua_out

    !! Work variables
    integer :: m


    ! Initialise - read VMR from table for T and p grid
    if (first_call .eqv. .True.) then
      call CE_interpolate_init(VMR_table)
      first_call = .False.
    end if

    lP_in = log10(P_in)
    lT_in = log10(T_in)

    ! Find upper and lower T and P triplet indexes
    call locate(P_t, nP, P_in, i_p2)
    i_p1 = i_p2 - 1
    i_p3 = i_p2 + 1

    if (i_p1 <= 0) then
      i_p1 = 1
      i_p2 = 2
      i_p3 = 3
    else if (i_p3 > nP) then
      i_p1 = nP - 2
      i_p2 = nP - 1
      i_p3 = nP
    end if

    ! Check if input temperature is within table range
    if (T_in <= T_t(1)) then

      ! Perform Bezier interpolation at minimum table temperature
      lPa(1) = lP_t(i_p1)
      lPa(2) = lP_t(i_p2)
      lPa(3) = lP_t(i_p3)
      mua(1) = mu_t(i_p1,1)
      mua(2) = mu_t(i_p2,1)
      mua(3) = mu_t(i_p3,1)
      call Bezier_interp(lPa(:), mua(:), 3, lP_in, mu_out)

      do m = 1, m_size
        xa(1) = x_t(m,i_p1,1)
        xa(2) = x_t(m,i_p2,1)
        xa(3) = x_t(m,i_p3,1)
        call Bezier_interp(lPa(:), xa(:), 3, lP_in, x_out(m))
        x_out(m) = 10.0_dp**x_out(m)
      end do

    else if (T_in >= T_t(nT)) then

      ! Perform Bezier interpolation at maximum table temperature
      lPa(1) = lP_t(i_p1)
      lPa(2) = lP_t(i_p2)
      lPa(3) = lP_t(i_p3)
      mua(1) = mu_t(i_p1,nT)
      mua(2) = mu_t(i_p2,nT)
      mua(3) = mu_t(i_p3,nT)
      call Bezier_interp(lPa(:), mua(:), 3, lP_in, mu_out)

      do m = 1, m_size
        xa(1) = x_t(m,i_p1,nT)
        xa(2) = x_t(m,i_p2,nT)
        xa(3) = x_t(m,i_p3,nT)
        call Bezier_interp(lPa(:), xa(:), 3, lP_in, x_out(m))
        x_out(m) = 10.0_dp**x_out(m)
      end do

    else

      ! Perform 2D Bezier interpolation by performing interpolation 4 times

      ! Find temperature index triplet
      call locate(T_t, nT, T_in, i_t2)
      i_t1 = i_t2 - 1
      i_t3 = i_t2 + 1

      lPa(1) = lP_t(i_p1)
      lPa(2) = lP_t(i_p2)
      lPa(3) = lP_t(i_p3)

      mua(1) = mu_t(i_p1,i_t1)
      mua(2) = mu_t(i_p2,i_t1)
      mua(3) = mu_t(i_p3,i_t1)
      call Bezier_interp(lPa(:), mua(:), 3, lP_in, mua_out(1)) ! Result at T1, P_in
      mua(1) = mu_t(i_p1,i_t2)
      mua(2) = mu_t(i_p2,i_t2)
      mua(3) = mu_t(i_p3,i_t2)
      call Bezier_interp(lPa(:), mua(:), 3, lP_in, mua_out(2)) ! Result at T2, P_in
      mua(1) = mu_t(i_p1,i_t3)
      mua(2) = mu_t(i_p2,i_t3)
      mua(3) = mu_t(i_p3,i_t3)
      call Bezier_interp(lPa(:), mua(:), 3, lP_in, mua_out(3)) ! Result at T3, P_in
      lTa(1) = lT_t(i_t1)
      lTa(2) = lT_t(i_t2)
      lTa(3) = lT_t(i_t3)
      call Bezier_interp(lTa(:), mua_out(:), 3, lT_in, mu_out) ! Result at T_in, P_in

      do m = 1, m_size
        xa(1) = x_t(m,i_p1,i_t1)
        xa(2) = x_t(m,i_p2,i_t1)
        xa(3) = x_t(m,i_p3,i_t1)
        call Bezier_interp(lPa(:), xa(:), 3, lP_in, xa_out(1)) ! Result at T1, P_in
        xa(1) = x_t(m,i_p1,i_t2)
        xa(2) = x_t(m,i_p2,i_t2)
        xa(3) = x_t(m,i_p3,i_t2)
        call Bezier_interp(lPa(:), xa(:), 3, lP_in, xa_out(2)) ! Result at T2, P_in
        xa(1) = x_t(m,i_p1,i_t3)
        xa(2) = x_t(m,i_p2,i_t3)
        xa(3) = x_t(m,i_p3,i_t3)
        call Bezier_interp(lPa(:), xa(:), 3, lP_in, xa_out(3)) ! Result at T3, P_in
        lTa(1) = lT_t(i_t1)
        lTa(2) = lT_t(i_t2)
        lTa(3) = lT_t(i_t3)
        call Bezier_interp(lTa(:), xa_out(:), 3, lT_in, x_out(m)) ! Result at T_in, P_in
        x_out(m) = 10.0_dp**x_out(m)
      end do

    end if

  end subroutine CE_interpolate_Bezier

  subroutine CE_interpolate_init(VMR_table)
    implicit none

    integer :: i, j
    character(len=50), intent(in) :: VMR_table

    !! Read T and P grid from file + mu and VMR data

    open(newunit=u, file=trim(VMR_table),action='read',form='formatted',status='old')

    read(u,*) nT, nP, npoint, nmol

    !print*, np, nmol

    allocate(m_name(nmol))
    read(u,*) (m_name(i),i=1,nmol)

    allocate(T_t(nT))
    read(u,*) (T_t(i),i=1,nT)

    allocate(P_t(nP))
    read(u,*) (P_t(i),i=1,nP)


    allocate(mu_t(nP,nT), x_t(nmol,nP,nT))
    do j = nP, 1,-1
      do i = nT, 1, -1
        read(u,*) mu_t(j,i), x_t(1:nmol,j,i)
      end do
    end do

    ! Log10 arrays of T-p grid
    allocate(lT_t(nT),lP_t(nP))
    lT_t = log10(T_t)
    lP_t = log10(P_t)

  end subroutine CE_interpolate_init

  subroutine locate(arr, n, var, idx)
    implicit none

    integer, intent(in) :: n
    integer, intent(out) :: idx
    real(dp), dimension(n), intent(in) :: arr
    real(dp), intent(in) ::  var
    integer :: jl, jm, ju

    ! Search an array using bi-section (numerical methods)
    ! Then return array index that is lower than var in arr

    jl = 0
    ju = n+1
    do while (ju-jl > 1)
      jm = (ju+jl)/2
      if ((arr(n) > arr(1)).eqv.(var > arr(jm))) then
        jl=jm
      else
        ju=jm
      end if
    end do

    idx = jl

  end subroutine locate

  subroutine linear_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    norm = 1.0_dp / (x2 - x1)

    yval = (y1 * (x2 - xval) + y2 * (xval - x1)) * norm

  end subroutine linear_interp

  pure subroutine bilinear_interp(xval, yval, x1, x2, y1, y2, a11, a21, a12, a22, aval)
    implicit none

    real(dp), intent(in) :: xval, yval, x1, x2, y1, y2, a11, a21, a12, a22
    real(dp), intent(out) :: aval
    real(dp) :: norm

    norm = 1.0_dp / (x2 - x1) / (y2 - y1)

    aval = a11 * (x2 - xval) * (y2 - yval) * norm &
      & + a21 * (xval - x1) * (y2 - yval) * norm &
      & + a12 * (x2 - xval) * (yval - y1) * norm &
      & + a22 * (xval - x1) * (yval - y1) * norm

  end subroutine bilinear_interp

  ! Perform Bezier interpolation
  subroutine Bezier_interp(xi, yi, ni, x, y)
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
      wlim = 1.0_dp + 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) - dx/2.0_dp * (w*dy/dx + (1.0_dp - w)*dy1/dx1)
      t = (x - xi(1))/dx
      y = (1.0_dp - t)**2 * yi(1) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(2)
    else ! (x > xi(2) and x < xi(3)) then
      ! right hand side interpolation
      !print*,'right'
      w = dx/(dx + dx1)
      wlim = 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp + 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) + dx1/2.0_dp * (w*dy1/dx1 + (1.0_dp - w)*dy/dx)
      t = (x - xi(2))/(dx1)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine Bezier_interp

end module CE_mod
