!!!
! Elspeth KH Lee - May 2021
! A simple program that emulates a column inside the Exo-FMS GCM
! This is useful for testing different RT solutions
! This version is for semi-grey radiation only
!
! Input parameters are set via the namelist (FMS_RC.nml) - hopefully paramaters are
! somewhat self explanatory. See the github readme for more information.
! Not guarenteed to be bug free
! NOTE: Indexing starts at 1, where 1 is the Top Of Atmosphere level or layer
!!!

program Exo_FMS_RC
  use, intrinsic :: iso_fortran_env

  use sw_direct_mod, only : sw_direct
  use sw_adding_mod, only : sw_adding
  use sw_SDA_mod, only : sw_SDA
  use sw_Toon_mod, only : sw_Toon
  use sw_SH2_mod, only : sw_SH2
  use sw_SH4_mod, only : sw_SH4
  use sw_disort_ts_mod, only : sw_disort_ts

  use lw_AA_E_mod, only : lw_AA_E
  use lw_AA_L_mod, only : lw_AA_L
  use lw_sc_linear_mod, only : lw_sc_linear
  use lw_VIM_mod, only : lw_VIM
  use lw_Toon_mod, only : lw_Toon
  use lw_disort_ts_mod, only : lw_disort_ts

  use ce_interp_mod, only : interp_ce_table
  use ce_Burrows_mod, only : analytic_Burrows

  use ck_opacity_mod, only : ck_opacity

  use IC_mod, only : IC_profile

  use dry_conv_adj_mod, only : Ray_dry_adj
  use MLT_mod, only : MLT

  use ieee_arithmetic
  implicit none

  ! Precision variable
  integer, parameter :: dp = REAL64

  ! Constants
  real(dp), parameter :: sb = 5.670374419e-8_dp
  real(dp), parameter :: R_gas = 8.31446261815324_dp
  real(dp), parameter :: Rsun = 6.95700e8_dp
  real(dp), parameter :: au = 1.495978707e11_dp

  integer :: n_ck, n_cia, n_ray
  character(len=200) ::  stellarf_sh, wl_sh

  integer :: n, i, k, u, j, b, inan
  integer :: nb, ng, nwl
  integer :: nstep, nlay, nlev
  integer :: table_num
  real(dp) :: t_step, t_tot
  real(dp) :: mu_z, Tirr, Tint, pref, pu, met, fl
  real(dp), allocatable, dimension(:) :: Tl, pl, pe, dpe, Finc
  real(dp), allocatable, dimension(:,:,:) :: k_l
  real(dp), allocatable, dimension(:,:,:) :: tau_e
  real(dp), allocatable, dimension(:,:,:) :: ssa, gg
  real(dp), allocatable, dimension(:) :: dT_rad, dT_conv, net_F
  real(dp), allocatable, dimension(:) :: gw
  real(dp), allocatable, dimension(:) :: wl_e, wn_e
  real(dp) :: olr, asr

  integer :: nsp
  character(len=200) :: VMR_tab_sh
  real(dp), allocatable, dimension(:,:) :: VMR
  real(dp), allocatable, dimension(:) :: mu
  character(len=10), allocatable, dimension(:) :: sp_list

  real(dp) :: grav, k_IR, k_V, cp_air, Rd_air, kappa_air
  real(dp) :: Rs, sm_ax

  real(dp), allocatable, dimension(:) :: cp_bar, kappa_bar, Rd_bar

  logical :: zcorr
  integer :: zcorr_meth
  real(dp) :: radius
  real(dp), allocatable, dimension(:) :: mu_z_eff, alt, alp

  real(dp), allocatable, dimension(:) :: a_surf
  real(dp), allocatable, dimension(:) :: sw_up, sw_down, sw_net 
  real(dp), allocatable, dimension(:) :: lw_up, lw_down, lw_net 

  real(dp), allocatable, dimension(:) :: Kzz 

  integer :: iIC
  logical :: corr
  real(dp) :: prc

  integer :: ua, ub, uu
  character(len=200) :: a_sh, b_sh
  real(dp), allocatable, dimension(:) :: a_hs, b_hs

  real(dp) :: start, finish

  character(len=200) :: sw_scheme, lw_scheme, opac_scheme, adj_scheme, CE_scheme

  integer :: u_nml

  namelist /FMS_RC_nml/ sw_scheme, lw_scheme, opac_scheme, adj_scheme, CE_scheme, nlay, a_sh, b_sh, pref, &
          & t_step, nstep, Rd_air, cp_air, grav, mu_z, Tirr, Tint, k_V, k_IR, fl, met, &
          & iIC, corr, table_num, nb, ng, nsp, wl_sh,  stellarf_sh, n_ck, &
          & n_cia, n_ray, Rs, sm_ax, zcorr, zcorr_meth, radius

  namelist /sp_nml/ sp_list, VMR_tab_sh

  namelist /gw_nml/ gw

  !! Read input variables from namelist
  open(newunit=u_nml, file='FMS_RC.nml', status='old', action='read')
  read(u_nml, nml=FMS_RC_nml)

  !! Read the species list and VMR table path
  allocate(sp_list(nsp))
  read(u_nml, nml=sp_nml)

  !! Read the Guassian ordinance weights for the corr-k opacity
  allocate(gw(ng))
  read(u_nml, nml=gw_nml)
  close(u_nml)

  !! Read wavelength and stellar flux in each band file
  allocate(Finc(nb))

  if (Tirr <= 0.0_dp) then
    Finc(:) = 0.0_dp
  else
    open(newunit=u,file=trim(stellarf_sh),status='old',action='read')
    do b = 1, nb
      read(u,*) Finc(b)       ! Read the stellar flux in each band
      !print*, b, Finc(b)
    end do
    close(u)

    Finc(:) = ((Rs * Rsun)/(sm_ax * au))**2 * Finc(:)
  end if

  open(newunit=u,file=trim(wl_sh),status='old',action='read')
  read(u,*) nwl
  allocate(wl_e(nwl), wn_e(nwl))
  if (nwl <= nb) then
    print*, 'Error, number of bands not compatable with number of wavelengths in file'
    stop
  end if
  do b = 1, nwl
    read(u,*) wl_e(b)
    wn_e(b) = 1.0_dp/(wl_e(b) * 1.0e-4_dp)
  end do
  close(u)

  !! Number of layer edges (levels)
  nlev = nlay + 1

  !! Read in hybrid sigma grid values
  open(newunit=ua,file=trim(a_sh), action='read', status='old')
  open(newunit=ub,file=trim(b_sh), action='read', status='old')
  allocate(a_hs(nlev),b_hs(nlev))
  do k = 1, nlev
    read(ua,*) a_hs(k)
    read(ub,*) b_hs(k)
  end do
  close(ua); close(ub)

  !! Contruct pressure array [pa] at the levels using the hybrid sigma formula
  ! Reference surface pressure [pa] is pref
  allocate(pe(nlev))
  do k = 1, nlev
    pe(k) = a_hs(k) + b_hs(k)*pref
  end do
  pu = pe(1)

  !! Pressure at the layers
  allocate(pl(nlay),dpe(nlay))
  do k = 1, nlay
    dpe(k) = pe(k+1) - pe(k)
    pl(k) = dpe(k) / log(pe(k+1)/pe(k))
  end do

  !! Allocate other arrays we need
  allocate(Tl(nlay), dT_rad(nlay), dT_conv(nlay), net_F(nlev))
  allocate(tau_e(ng,nb,nlev), k_l(ng,nb,nlay))
  allocate(ssa(ng,nb,nlay), gg(ng,nb,nlay))
  allocate(VMR(nsp,nlay), mu(nlay))
  allocate(alt(nlev), mu_z_eff(nlev), alp(nlev))
  allocate(sw_up(nlev), sw_down(nlev), sw_net(nlev))
  allocate(lw_up(nlev), lw_down(nlev), lw_net(nlev))
  allocate(a_surf(nb))
  allocate(Kzz(nlay))
  allocate(cp_bar(nlay), kappa_bar(nlay), Rd_bar(nlay))

  !! Set a_surf to zero for now
  a_surf(:) = 0.0_dp

  !! Make thermo variables constabt for now
  cp_bar(:) = cp_air
  Rd_bar(:) = Rd_air

  !! Calculate the adiabatic coefficent
  kappa_air = Rd_air/cp_air   ! kappa = Rd/cp
  kappa_bar(:) = Rd_bar(:)/cp_bar(:)

  print*, 'Tint ', 'Tirr ', 'pref ', 'pu ', 'mu_z ', 'grav '
  print*, Tint, Tirr, pref/1e5_dp, pu/1e5_dp, mu_z, grav
  print*, '--------------'

  ! Semi-grey atmosphere values (here they are not used, but just need to be passed to IC routine)
  k_l(1,1,:) = k_V
  k_l(1,1,:) = k_IR

  fl = 1.0_dp

  !! Initial condition T-p profile - see the routine for options
  call IC_profile(iIC,corr,nlay,pref,pl,k_l(1,1,:),k_l(1,1,:),Tint,mu_z,Tirr,grav,fl,Tl,prc,table_num,met)

  !! Print initial T-p profile
  do i = 1, nlay
    print*, i, pl(i)/1e5_dp, Tl(i)
  end do
  print*, '--------------'

  !! Write out initial conditions
  open(newunit=u,file='FMS_RC_ic.out',action='readwrite')
  do i = 1, nlay
    write(u,*) i, pl(i), Tl(i)
  end do
  close(u)

  !! Time stepping loop
  print*, 'Start timestepping'

  t_tot = 0.0_dp
  inan = 0

  !! cpu timer start
  call cpu_time(start)

  do n = 1, nstep

    net_F(:) = 0.0_dp
    dT_conv(:) = 0.0_dp

    select case(CE_scheme)
    case('interp')
      do k = 1, nlay
        call interp_ce_table(nsp, Tl(k), pl(k), VMR(:,k), mu(k), VMR_tab_sh)
        !print*, Tl(k), pl(k)/1e5_dp, mu(k), VMR(:,k)
      end do
    case('Burrows')
      do k = 1, nlay
        call analytic_Burrows(nsp, Tl(k), pl(k), VMR(:,k))
        mu(k) = R_gas/Rd_air * 1000.0_dp
      end do
    case('none')
      do k = 1, nlay
        mu(k) = R_gas/Rd_air * 1000.0_dp
      end do
    case default
      print*, 'Invalid CE_scheme: ', trim(CE_scheme)
      stop
    end select

    select case(opac_scheme)
    case('ck')

       ! Calculate the opacity structure of the column
       call ck_opacity(nlay, n_ck, n_CIA, n_Ray, nb, ng, wl_e, grav, Tl(:), pl(:), pe(:), mu(:), &
         & nsp, sp_list(:), VMR(:,:), k_l(:,:,:), ssa(:,:,:), gg(:,:,:))

      ! Include optical depth component from 0 pressure, assuming constant T and p at boundary
      tau_e(:,:,1) = (k_l(:,:,1) * pe(1)) / grav
      do k = 1, nlay
        tau_e(:,:,k+1) = tau_e(:,:,k) + (k_l(:,:,k) * dpe(k)) / grav
        !print*, k, k_l(1,1,k), tau_e(1,1,k+1), k_l(8,1,k), tau_e(8,1,k+1)
      end do

    case default
      print*, 'Invalid opac_scheme: ', trim(opac_scheme)
      stop
    end select

    !! Zenith angle geometric correction
    if ((zcorr .eqv. .True.) .and. (mu_z > 0.0_dp)) then
      ! First calculate the altitude at each level from the hypsometric equation
      ! Assume constant gravity for simplicity
      alt(nlev) = 0.0_dp
      do k = nlay, 1, -1
        alt(k) = alt(k+1) + (Rd_air*Tl(k))/grav * log(pe(k+1)/pe(k))
      end do

      select case(zcorr_meth)
      case(1)
        ! Basic geometric correction following Li & Shibata (2006) Eq. (2)
        mu_z_eff(:) = sqrt(1.0_dp - (radius/(radius + alt(:)))**2 * (1.0_dp - mu_z**2))
      case(2)
        ! Spherical layer correction following Li & Shibata (2006) Eq.(10)
        alp(nlev) = (alt(nlay) -  alt(nlev))/radius
        do k = nlay,1,-1
           alp(k) = (alt(k) -  alt(k+1))/(radius + alt(k))
        end do
        mu_z_eff(:) = (sqrt(1.0_dp - (radius/(radius + alt(:)))**2 * (1.0_dp - mu_z**2)) + &
          & sqrt((1.0_dp + alp(:))**2 - (radius/(radius + alt(:)))**2 * (1.0_dp - mu_z**2))) / &
          & (2.0_dp + alp(:))
      case default
        print*, 'Invalid zcorr_meth ', zcorr_meth
        stop
      end select
    else
      ! No correction, use single zenith angle
      mu_z_eff(:) = mu_z
    end if

    !! Shortwave radiative transfer step
    select case(sw_scheme)
    case('sw_direct')
      ! Direct beam only method with no scattering
      call sw_direct(nlev, nb, ng, gw, tau_e, mu_z_eff, Finc, a_surf, sw_up, sw_down, sw_net, asr)
    case('sw_adding')
      ! Approximate adding method with approximate scattering
      call sw_adding(nlay, nlev, nb, ng, gw, tau_e, mu_z_eff, Finc, ssa, gg, a_surf, sw_up, sw_down, sw_net, asr)
    case('sw_SDA')
      ! Spherical harmonic doubling (SDA) adding four stream method with approximate scattering
      call sw_SDA(nlay, nlev, nb, ng, gw, tau_e, mu_z_eff, Finc, ssa, gg, a_surf, sw_up, sw_down, sw_net, asr)
    case('sw_Toon')
      ! Toon89 shortwave method with multiple scattering
      call sw_Toon(nlay, nlev, nb, ng, gw, tau_e, mu_z_eff, Finc, ssa, gg, a_surf, sw_up, sw_down, sw_net, asr)
    case('sw_SH2')
      ! Spherical harmonic two stream method with multiple scattering
      call sw_SH2(nlay, nlev, nb, ng, gw, tau_e, mu_z_eff, Finc, ssa, gg, a_surf, sw_up, sw_down, sw_net, asr)
    case('sw_SH4')
      ! Spherical harmonic four stream method with multiple scattering
      call sw_SH4(nlay, nlev, nb, ng, gw, tau_e, mu_z_eff, Finc, ssa, gg, a_surf, sw_up, sw_down, sw_net, asr)
    case('sw_disort_ts')
      ! Two stream disort method with multiple scattering
      call sw_disort_ts(nlay, nlev, nb, ng, gw, tau_e, mu_z_eff, Finc, ssa, gg, sw_up, sw_down, sw_net, asr)
    case('none')
    case default
      print*, 'Invalid sw_scheme: ', trim(sw_scheme)
      stop
    end select

    !! Longwave radiative transfer step
    select case(lw_scheme)
    case('lw_AA_E')
      ! Absorption approximation exponential in tau (AA_E) with approximate scattering
      call lw_AA_E(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr) 
    case('lw_AA_L')
      ! Absorption approximation linear in tau (AA_L) with approximate scattering
      call lw_AA_L(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)     
    case('lw_sc_linear')
      ! Short characteristics (sc) with linear interpolants with no scattering
      call lw_sc_linear(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_VIM')
      ! Variational Iteration Method (VIM) with approximate scattering
      call lw_VIM(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_Toon')
      ! Toon89 longwave method with multiple scattering
      call lw_Toon(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, a_surf, Tint, lw_up, lw_down, lw_net, olr)
    case('lw_disort_ts')
      ! Two stream disort method with multiple scattering
      call lw_disort_ts(nlay, nlev, nb, ng, gw, wn_e, Tl, pl, pe, tau_e, ssa, gg, Tint, lw_up, lw_down, lw_net, olr) 
    case('none')
    case default
      print*, 'Invalid lw_scheme: ', trim(lw_scheme)
      stop
    end select


    !! Calculate the temperature tendency due to radiation
    net_F(:) = lw_net(:) + sw_net(:) ! Net flux into/out of level
    do i = 1, nlay
      dT_rad(i) = (grav/cp_air) * (net_F(i+1)-net_F(i))/(dpe(i))
    end do

    !! Convective adjustment scheme
    select case(adj_scheme)
    case('Ray_dry')
      ! Dry convective adjustment following Ray Pierrehumbert's python script
      call Ray_dry_adj(nlay, nlev, t_step, kappa_air, Tl, pl, pe, dT_conv)
      Kzz(:) = 1e1_dp
    case('MLT')
      ! Use mixing length theory (MLT) to time dependently adjust the adiabat and estimate Kzz
      call MLT(nlay, nlev, t_step, Tl, pl, pe, Rd_bar, cp_bar, kappa_bar, &
         & grav, dT_conv, Kzz)      
    case('none')
    case default
      print*, 'Invalid adj_scheme: ', trim(adj_scheme)
      stop
    end select

    !! Forward march the temperature change in each layer from convection and radiation
    Tl(:) = Tl(:) + t_step * (dT_conv(:) + dT_rad(:))

    !! Check for NaN's in the temperature and exit the simulation if detected
    do i = 1, nlay
      if (ieee_is_nan(Tl(i)) .eqv. .True.) then
        do j = 1, nlay
          print*, j, Tl(j), net_F(j), dT_rad(j), dT_conv(j)
        end do
        print*, nlev, net_F(nlev)
        inan = 1
        exit
      end if
    end do
    if (inan == 1) then
      exit
    end if

    !! Increase the total time simulated
    t_tot = t_tot + t_step

  end do

  !! cpu timer end
  call cpu_time(finish)

  !! Output the results
  print*, 'sec: ', 'hours: ', 'days: '
  print*, t_tot, t_tot/60.0_dp/60.0_dp,t_tot/60.0_dp/60.0_dp/24.0_dp

  print*, 'For profile properties: '
  print*, Tint, Tirr, pref, mu_z

  print*, 'OLR [W m-2], Teff:'
  print*, olr, (olr/sb)**(0.25_dp)

  print*, 'ASR [W m-2], Tinc:'
  print*, asr, (asr/sb)**(0.25_dp)

  print*, 'Internal T [W m-2], Tint'
  print*, sb * Tint**4, Tint

  print*, 'Outputting results: '
  open(newunit=uu,file='FMS_RC_pp.out', action='readwrite')
  do i = 1, nlay
    write(uu,*) i, pl(i), Tl(i), dT_rad(i), dT_conv(i), Kzz(i)
  end do
  close(uu)

  open(newunit=uu,file='FMS_RC_flx.out', action='readwrite')
  do i = 1, nlev
    write(uu,*) i, pe(i), sw_up(i), sw_down(i), sw_net(i), lw_up(i), lw_down(i), lw_net(i)
  end do
  close(uu)

  print*, n, 'steps took: '
  print '("Time = ",f8.3," seconds.")', finish-start


end program Exo_FMS_RC
