!! This module

module IC_mod
  use, intrinsic :: iso_fortran_env
  use k_Rosseland_mod, only : k_Ross_Freedman, k_Ross_Valencia, gam_Parmentier, Bond_Parmentier
  implicit none

  !! Precision variables
  integer, parameter :: dp = REAL64

  public :: IC_profile
  private :: adiabat_correction, Parmentier_IC, Guillot_IC, Tdeep_IC, Iso_IC, Mayne_2014_IC

contains

  subroutine IC_profile(iIC,corr,nlay,p0,pl,k_V,k_IR,Tint,mu,Tirr,grav,fl,Tl,prc,table_num,met)
    implicit none

    !! Input flags
    integer, intent(in) :: iIC
    logical, intent(in) :: corr

    !! Input quantities
    integer, intent(in) :: nlay, table_num
    real(dp), intent(in) :: p0, Tint, mu, Tirr, grav, fl, met
    real(dp), dimension(nlay), intent(in) :: pl
    real(dp), dimension(nlay), intent(in) :: k_V, k_IR

    !! Output quantities
    real(dp), dimension(nlay), intent(out) :: Tl
    real(dp), intent(out) :: prc

    select case(iIC)

    case(1)
      ! Isothermal IC
      call Iso_IC(nlay,Tint,Tl)
    case(2)
      ! Deep temperature IC following Guillot (2010)
      call Tdeep_IC(nlay,k_V(1),k_IR(1),Tirr,Tl)
    case(3)
      ! Guillot (2010) analytic RE T-p profile
      call Guillot_IC(nlay,p0,pl,k_V(1),k_IR(1),Tint,mu,Tirr,grav,fl,Tl)
    case(4)
      ! Parmentier et al. (2014, 2015) IC picket fence IC
      ! check the table_num and metallicity parameters for different options
      call Parmentier_IC(nlay,pl,Tint,mu,Tirr,grav,Tl,table_num,met)
    case(5)
      ! Mayne et al. (2014) IC - check the iprof parameter in the subroutine
      ! for then dayside/nightside switch
      call Mayne_2014_IC(nlay,pl,Tl)
    case(6)
      call CAMEMBERT_IC(nlay,pl,Tl)
    case default
      print*, 'Invalid IC integer in IC_mod, stopping'
      stop

    end select

    !! Perform adiabatic correction according to Parmentier et al. (2015)
    if (corr .eqv. .True.) then
      call adiabat_correction(nlay,Tl,pl,prc)
    else
      prc = p0
    end if

  end subroutine IC_profile

  subroutine CAMEMBERT_IC(nlay,pl,Tl)
    implicit none

    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: pl

    real(dp), dimension(nlay), intent(out) :: Tl

    real(dp), dimension(:), allocatable :: p_C, T_C

    integer :: i, nTP, u

    integer :: iT, iT1, ip, ip1

    open(newunit=u,file='CAMEMBERT_GJ1214b_IC_H2He.dat',action='read')

    do i = 1, 6
      read(u,*)
    end do

    nTP = 87

    allocate(p_C(nTP), T_C(nTP))

    do i = 1, nTP
      read(u,*) p_C(i), T_C(i)
    end do

    do i = 1, nlay

      call locate(p_C,pl(i),ip)
      ip1 = ip + 1

      if (ip == 0) then
        Tl(i) = T_C(1)
      else if (ip == nTP) then
        Tl(i) = T_C(nTP)
      else
        call linear_log_interp(pl(i), p_C(ip), p_C(ip1), T_C(ip), T_C(ip1), Tl(i))
      end if

    end do

  end subroutine CAMEMBERT_IC

    subroutine locate(arr, var, idx)
    implicit none

    integer, intent(out) :: idx
    real(kind=dp), dimension(:), intent(in) :: arr
    real(kind=dp),intent(in) ::  var
    integer :: jl, jm, ju

    ! Search an array using bi-section/binary search (numerical methods)
    ! Then return array index that is lower than var in arr

    jl = 0
    ju = size(arr)+1
    do while (ju-jl > 1)
      jm = (ju+jl)/2
      if (var > arr(jm)) then
        jl=jm
      else
        ju=jm
      end if
    end do

    idx = jl

  end subroutine locate

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(kind=dp), intent(in) :: xval, y1, y2, x1, x2
    real(kind=dp) :: lxval, ly1, ly2, lx1, lx2
    real(kind=dp), intent(out) :: yval
    real(kind=dp) :: norm

    lxval = log10(xval)
    lx1 = log10(x1); lx2 = log10(x2)
    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / (lx2 - lx1)

    yval = 10.0_dp**((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)

  end subroutine linear_log_interp

  subroutine Iso_IC(nlay,Tint,Tl)
    implicit none

    !! Input quantities
    integer, intent(in) :: nlay
    real(dp), intent(in) :: Tint

    !! Output quantities
    real(dp), dimension(nlay), intent(out) :: Tl

    Tl(:) = Tint

  end subroutine Iso_IC

  subroutine Tdeep_IC(nlay,k_V,k_IR,Tirr,Tl)
    implicit none

    !! Input quantities
    integer, intent(in) :: nlay
    real(dp), intent(in) :: k_V, k_IR, Tirr

    !! Output quantities
    real(dp), dimension(nlay), intent(out) :: Tl

    !! Work variables
    real(dp) :: gam

    gam = k_V/k_IR

    Tl(:) = Tirr * (3.0_dp/4.0_dp * (1.0_dp/gam + 2.0_dp/3.0_dp))**(1.0_dp/4.0_dp)

  end subroutine Tdeep_IC

  subroutine Guillot_IC(nlay,p0,pl,k_V,k_IR,Tint,mu,Tirr,grav,fl,Tl)
    implicit none

    !! Input quantities
    integer, intent(in) :: nlay
    real(dp), intent(in) :: p0, k_V, k_IR, Tint, mu, Tirr, grav, fl
    real(dp), dimension(nlay), intent(in) :: pl

    !! Output quantities
    real(dp), dimension(nlay), intent(out) :: Tl

    !! Work variables
    real(dp), dimension(nlay) :: tau_IRl
    real(dp) :: gam, tau0

    gam = k_V/k_IR
    tau0 = k_IR/grav * p0 / fl

    tau_IRl(:) = fl * (pl(:)/p0 * tau0) + (1.0_dp - fl) * ((pl(:)/p0)**2 * tau0)

    Tl(:) = ((3.0_dp/4.0_dp) * Tint**4 * (tau_IRl(:) + 2.0_dp/3.0_dp))
    Tl(:) = Tl(:) + (mu * 3.0_dp * Tirr**4)/4.0_dp *  &
         & (2.0_dp/3.0_dp + mu/gam + ((gam/(3.0_dp*mu)) - mu/gam) * exp(-gam*tau_IRl(:)/mu))
    Tl(:) = Tl(:)**(1.0_dp/4.0_dp)

  end subroutine Guillot_IC

  !! This subroutine follows Parmentier & Guillot (2014, 2015) non-grey picket fence scheme
  subroutine Parmentier_IC(nlay,pl,Tint,mu,Tirr,grav,Tl,table_num,met)
    implicit none

    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: pl
    real(dp), intent(in) :: Tint, mu, Tirr, grav
    integer, intent(in) :: table_num
    real(dp), intent(in) :: met


    real(dp), dimension(nlay), intent(out) :: Tl

    real(dp) :: Teff0, Teff, Tmu, Bond, Tskin
    real(dp), dimension(3) :: gam_V, Beta_V
    real(dp), dimension(2) :: Beta
    real(dp) :: gam_1, gam_2, gam_P, tau_lim

    integer :: i, j
    real(dp) :: a0, a1, b0, A, B, At1, At2
    real(dp), dimension(3) :: a2, a3, b1, b2, b3, Av1, Av2
    real(dp), dimension(3) :: C, D, E
    real(dp), dimension(nlay+1) :: tau
    real(dp), dimension(nlay) :: kRoss

    !! Effective temperature parameter
    Tmu = (mu * Tirr**4)**(1.0_dp/4.0_dp)

    !! Find Bond albedo of planet - Bond albedo is given by mu = 1/sqrt(3)
    Teff0 = (Tint**4 + (1.0_dp/sqrt(3.0_dp)) * Tirr**4)**(1.0_dp/4.0_dp)
    call Bond_Parmentier(Teff0, grav, Bond)

    Teff = (Tint**4 + (1.0_dp - Bond) * mu * Tirr**4)**(1.0_dp/4.0_dp)

    !! Find the V band gamma, beta and IR gamma and beta ratios for this profile
    ! Passed mu, so make lat = acos(mu) and lon = 0
    call gam_Parmentier(Teff, table_num, gam_V, Beta_V, Beta, gam_1, gam_2, gam_P, tau_lim)

    gam_V(:) = gam_V(:) / mu

    !! Hard work starts here - first calculate all the required coefficents
    At1 = gam_1**2*log(1.0_dp + 1.0_dp/(tau_lim*gam_1))
    At2 = gam_2**2*log(1.0_dp + 1.0_dp/(tau_lim*gam_2))
    Av1(:) = gam_1**2*log(1.0_dp + gam_V(:)/gam_1)
    Av2(:) = gam_2**2*log(1.0_dp + gam_V(:)/gam_2)

    a0 = 1.0_dp/gam_1 + 1.0_dp/gam_2

    a1 = -1.0_dp/(3.0_dp * tau_lim**2) * (gam_P/(1.0_dp-gam_P) * (gam_1 + gam_2 - 2.0_dp)/(gam_1 + gam_2) &
    &  + (gam_1 + gam_2)*tau_lim - (At1 + At2)*tau_lim**2)

    a2(:) = tau_lim**2/(gam_P*gam_V(:)**2) &
    &  * ((3.0_dp*gam_1**2-gam_V(:)**2)*(3.0_dp*gam_2**2-gam_V(:)**2)*(gam_1+gam_2) &
    & - 3.0_dp*gam_V(:)*(6.0_dp*gam_1**2*gam_2**2-gam_V(:)**2*(gam_1**2+gam_2**2))) &
    & / (1.0_dp-gam_V(:)**2 * tau_lim**2)

    a3(:) = -tau_lim**2*(3.0_dp*gam_1**2-gam_V(:)**2)*(3.0_dp*gam_2**2-gam_V(:)**2)*(Av2(:)+Av1(:)) &
     &/(gam_P*gam_V(:)**3*(1.0_dp-gam_V(:)**2*tau_lim**2))

    b0 = 1.0_dp/(gam_1*gam_2/(gam_1-gam_2)*(At1-At2)/3.0_dp-(gam_1*gam_2)**2/sqrt(3.0_dp*gam_P)-(gam_1*gam_2)**3 &
    & / ((1.0_dp-gam_1)*(1.0_dp-gam_2)*(gam_1+gam_2)))

    b1(:) = gam_1*gam_2*(3.0_dp*gam_1**2-gam_V(:)**2)*(3.0_dp*gam_2**2-gam_V(:)**2)*tau_lim**2 &
    & / (gam_P*gam_V(:)**2*(gam_V(:)**2*tau_lim**2-1.0_dp))

    b2(:) = 3.0_dp*(gam_1+gam_2)*gam_V(:)**3/((3.0_dp*gam_1**2-gam_V(:)**2)*(3.0_dp*gam_2**2-gam_V(:)**2))

    b3(:) = (Av2(:)-Av1(:))/(gam_V(:)*(gam_1-gam_2))

    A = 1.0_dp/3.0_dp*(a0+a1*b0)
    B = -1.0_dp/3.0_dp*(gam_1*gam_2)**2/gam_P*b0
    C(:) = -1.0/3.0_dp*(b0*b1(:)*(1.0_dp+b2(:)+b3(:))*a1+a2(:)+a3(:))
    D(:) = 1.0/3.0_dp*(gam_1*gam_2)**2/gam_P*b0*b1(:)*(1.0_dp+b2(:)+b3(:))
    E(:) = (3.0_dp-(gam_V(:)/gam_1)**2)*(3.0_dp-(gam_V(:)/gam_2)**2)/(9.0_dp*gam_V(:)*((gam_V(:)*tau_lim)**2-1.0_dp))

    ! T-p structure calculation - we follow exactly V. Parmentier's method
    ! Estimate the skin temperature by setting tau = 0
    tau(1) = 0.0_dp
    Tskin = 3.0_dp*Tint**4/4.0_dp*(tau(1)+A+B*exp(-tau(1)/tau_lim)) + sum(3.0_dp*Beta_V(:) &
    & * Tmu**4/4.0_dp*(C(:)+D(:)*exp(-tau(1)/tau_lim)+E(:)*exp(-gam_V(:)*tau(1))))
    Tskin = Tskin**(1.0_dp/4.0_dp)
    ! Estimate the opacity TOA at the skin temperature - assume this is = first layer optacity
    call k_Ross_Freedman(Tskin, pl(1), met, kRoss(1))
    !call k_Ross_Valencia(Tskin, pe(1), met, kRoss(1))

    ! Recalculate the upmost tau with new kappa
    tau(1) = kRoss(1)/grav * pl(1)
    ! More accurate layer T at uppermost layer
    Tl(1) = 3.0_dp*Tint**4/4.0_dp*(tau(1)+A+B*exp(-tau(1)/tau_lim)) + sum(3.0_dp*Beta_V(:) &
    & * Tmu**4/4.0_dp*(C(:)+D(:)*exp(-tau(1)/tau_lim)+E(:)*exp(-gam_V(:)*tau(1))))
    Tl(1) = Tl(1)**(1.0_dp/4.0_dp)

    ! Now we can loop in optical depth space to find the T-p profile
    do i = 2, nlay
      ! Initial guess for layer
      call k_Ross_Freedman(Tl(i-1), sqrt(pl(i-1)*pl(i)), met, kRoss(i))
      !call k_Ross_Valencia(Tl(i-1), sqrt(pl(i-1)*pl(i)), met, kRoss(i))
      tau(i) = tau(i-1) + kRoss(i)/grav * (pl(i) - pl(i-1))
      Tl(i) = 3.0_dp*Tint**4/4.0_dp*(tau(i)+A+B*exp(-tau(i)/tau_lim)) + sum(3.0_dp*Beta_V(:) &
      & * Tmu**4/4.0_dp*(C(:)+D(:)*exp(-tau(i)/tau_lim)+E(:)*exp(-gam_V(:)*tau(i))))
      Tl(i) = Tl(i)**(1.0_dp/4.0_dp)
      ! Convergence loop
      do j = 1, 5
        call k_Ross_Freedman(sqrt(Tl(i-1)*Tl(i)), sqrt(pl(i-1)*pl(i)), met, kRoss(i))
        !call k_Ross_Valencia(sqrt(Tl(i-1)*Tl(i)), sqrt(pl(i-1)*pl(i)), met, kRoss(i))
        tau(i) = tau(i-1) + kRoss(i)/grav * (pl(i) - pl(i-1))
        Tl(i) = 3.0_dp*Tint**4/4.0_dp*(tau(i)+A+B*exp(-tau(i)/tau_lim)) + sum(3.0_dp*Beta_V(:) &
        & * Tmu**4/4.0_dp*(C(:)+D(:)*exp(-tau(i)/tau_lim)+E(:)*exp(-gam_V(:)*tau(i))))
        Tl(i) = Tl(i)**(1.0_dp/4.0_dp)

      end do

    end do

  end subroutine Parmentier_IC

  !! Subroutine that corrects for adiabatic region following Parmentier & Guillot (2015)
  subroutine adiabat_correction(nlay,Tl,pl,prc)
    implicit none

    !! Input quantities
    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) ::  pl

    !! Output quantities
    real(dp), dimension(nlay), intent(inout) :: Tl
    real(dp), intent(out) :: prc

    !! Work variables
    integer :: k, iRC, iRC1
    real(dp), dimension(nlay) :: gradrad, gradad

    do k = 1, nlay-1
      gradrad(k) = (log10(Tl(k))-log10(Tl(k+1)))/(log10(pl(k))-log10(pl(k+1)))
      gradad(k) = 0.32_dp - 0.10_dp*Tl(k)/3000.0_dp
      !print*, k, gradrad(k), gradad(k)
    end do
    gradrad(nlay) = 0.0_dp
    gradad(nlay) = 0.0_dp

    iRC = nlay-1
    iRC1 = nlay-1

    do k = nlay-1, 1, -1
      if (IRC1 <= k+1) then
        if (gradrad(k) > 0.7_dp*gradad(k)) then
          iRC1 = k
        endif
        if (gradrad(k) > 0.98_dp*gradad(k)) then
         iRC = k
         prc = pl(iRC)
        endif
      end if
    end do

    if (iRC < nlay) then
      do k = iRC, nlay-1
        gradad(k) = 0.32_dp - 0.10_dp*Tl(k)/3000.0_dp
        if (gradad(k) < 0.0_dp) then
          gradad(k) = 0.0_dp
        end if
        !if (pl(k) > 1.0_dp*1e6_dp) then
        Tl(k+1)=Tl(k)*(pl(k+1)/pl(k))**gradad(k)
      !end if
      end do
    end if

  end subroutine adiabat_correction

  subroutine Mayne_2014_IC(nlay,pl,Tl)
    implicit none

    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: pl
    real(dp), dimension(nlay), intent(out) :: Tl

    integer :: i, iprof
    real(dp) :: p_low, p_hi, pl_b, lpl_b
    real(dp) :: Tp_day, T_day, Tp_night, T_night

    iprof = 1
    p_low = 100.0_dp
    p_hi = 1e6_dp

    if (iprof == 1) then
      ! Nightside profile
      do i = 1, nlay
        if (pl(i) < p_low) then
          pl_b = p_low / 1.0e5_dp
        else if (pl(i) >= p_hi) then
          pl_b = p_hi / 1.0e5_dp
        else
          pl_b = pl(i) / 1.0e5_dp
        end if
        lpl_b = log10(pl_b)
        if (pl_b <= 10.0_dp) then
          Tp_night = 1388.2145_dp + 267.66586_dp*lpl_b - 215.53357_dp*lpl_b**2 + 61.814807_dp*lpl_b**3 + 135.68661_dp*lpl_b**4 &
          & + 2.0149044_dp*lpl_b**5 - 40.907246_dp*lpl_b**6 - 19.015628_dp*lpl_b**7 - 3.8771634_dp*lpl_b**8 - &
          & 0.38413901_dp*lpl_b**9 - 0.015089084_dp*lpl_b**10
        else
          Tp_night = 5529.7168_dp - 6869.6504_dp*lpl_b + 4142.7231_dp*lpl_b**2 - 936.23053_dp*lpl_b**3 + 87.120975_dp*lpl_b**4
        end if

        if (pl(i) >= p_hi) then
          T_night = Tp_night + 100.0_dp * (1.0_dp - exp(-(log10(pl(i)) - log10(p_hi))))
        else if (pl(i) < p_low) then
          T_night = max(Tp_night * exp(0.10_dp*(log10(pl(i)) - log10(p_low))), 250.0_dp)
        else
          T_night = Tp_night
        end if
        Tl(i) = T_night
      end do
    else if (iprof == 2) then
      ! dayside profile
      do i = 1, nlay
        if (pl(i) < p_low) then
          pl_b = p_low / 1.0e5_dp
        else if (pl(i) >= p_hi) then
          pl_b = p_hi / 1.0e5_dp
        else
          pl_b = pl(i) / 1.0e5_dp
        end if
        lpl_b = log10(pl_b)
        if (pl_b <= 10.0_dp) then
          Tp_day = 2149.9581_dp + 4.1395571_dp*lpl_b - 186.24851_dp*lpl_b**2 + 135.525_dp*lpl_b**3 + 106.20433_dp*lpl_b**4 &
          & - 35.851966_dp*lpl_b**5 - 50.022826_dp*lpl_b**6 - 18.462489_dp*lpl_b**7 - 3.3319965_dp*lpl_b**8 - &
          & 0.30295925_dp*lpl_b**9 - 0.011122316_dp*lpl_b**10
        else
          Tp_day = 5529.7168_dp - 6869.6504_dp*lpl_b + 4142.7231_dp*lpl_b**2 - 936.23053_dp*lpl_b**3 + 87.120975_dp*lpl_b**4
        end if

        if (pl(i) >= p_hi) then
          T_day = Tp_day - 120.0_dp * (1.0_dp - exp(-(log10(pl(i)) - log10(p_hi))))
        else if (pl(i) < p_low) then
          T_day = max(Tp_day * exp(0.015_dp*(log10(pl(i)) - log10(p_low))), 1000.0_dp)
        else
          T_day = Tp_day
        end if
        Tl(i) = T_day
      end do

    end if


  end subroutine Mayne_2014_IC


end module IC_mod

! program IC_mod_test
!   use IC_mod, only : IC_profile
!   use, intrinsic :: iso_fortran_env
!   implicit none
!
!   integer, parameter :: dp = REAL64
!
!   integer :: iIC, nlay, nlay1, k
!   logical :: corr
!
!   real(dp) :: p0, mu, Tint, Tirr, fl, grav, prc, Teff
!   real(dp), allocatable, dimension(:) :: pe, pl, Tl, k_V_l, k_IR_l
!
!   integer :: u, ua, ub
!   character(len=20) :: a_sh, b_sh
!   real(dp), allocatable, dimension(:) :: a, b
!
!   ! Number of layers and edges
!   nlay = 52
!   nlay1 = nlay + 1
!
!   !! Read in sigma hybrid grid values
!   a_sh = 'sig_hyb_HJ_53_a.txt'
!   b_sh = 'sig_hyb_HJ_53_b.txt'
!
!   open(newunit=ua,file=trim(a_sh),action='read')
!   open(newunit=ub,file=trim(b_sh),action='read')
!
!   allocate(a(nlay1),b(nlay1))
!   do k = 1, nlay1
!     read(ua,*) a(k)
!     read(ub,*) b(k)
!   end do
!
!   ! Contruct pressure array in pa
!   ! Surface pressure (pa)
!   p0 = 1e8_dp
!   allocate(pe(nlay1),pl(nlay))
!   do k = 1, nlay1
!     pe(k) = a(k) + b(k)*p0
!     !print*, pe(k)/1e5_dp
!   end do
!   ! Pressure layers
!   do k = 1, nlay
!     pl(k) = (pe(k) + pe(k+1)) / 2.0_dp
!   end do
!
!   allocate(Tl(nlay))
!   allocate(k_V_l(nlay),k_IR_l(nlay))
!
!   ! sw Zenith angle
!   mu = 1.0_dp / sqrt(3.0_dp)
!
!   Tirr = 3335.0_dp ! Irradiation temperature
!   Tint = 583.0_dp ! Internal temperature
!
!   Teff = (Tint**4 + mu*Tirr**4)**(1.0_dp/4.0_dp)
!
!   k_IR_l(:) = 1e-3_dp
!   k_V_l(:) = 6e-4_dp * sqrt(Tirr/2000.0_dp)
!
!   grav = 10.8_dp
!   fl = 1.0_dp
!   Tl = 0.0_dp
!
!   iIC = 4
!   corr = .True.
!
!   call IC_profile(iIC,corr,nlay,p0,pl,pe,k_V_l,k_IR_l,Tint,mu,Tirr,grav,fl,Tl,prc)
!
!   open(newunit=u,file='FMS_IC.out',action='readwrite')
!   do k = 1, nlay
!     write(u,*) k, pl(k), Tl(k), prc
!   end do
!   close(u)
!
!   print*, Tirr, Tint, Teff, grav, mu,  (mu*Tirr**4)**(1.0_dp/4.0_dp), prc/1e5
!
!
! end program IC_mod_test
