!!!
! Elspeth K.H. Lee (Jun 2024)
! A large module for calculating corr-k, cia and rayleigh opacities
!   - Nov 2026 : Flux-weighted equivalent extinction (EE) and the RORR rebinning fix, following
!                Exo-FMS_2 (exp/*/srcmods/ck_opacity_mod.f90)
!   - Nov 2026 : One shared monotonic Bezier interpolation (bezier_index, bezier_1D, bezier_2D,
!                Bezier_interp - identical to Exo-FMS_2) for every interpolated opacity source
!                with ck_interp = 2: k-tables in log10 T, log10 p; CIA form 3 in log10 T; CIA form 1
!                in wn and log10 T; free-free tables (form 2) in log10 wl and log10 theta; all on
!                log10 values, edge values outside the tables. The free-free tables now use their
!                edge values outside the table (they gave zero below T = 5040/3.6 = 1400 K).
!
! k-table mixing methods (ck_nml, choose one):
!   PM   : a premixed table (one table, already per gram)
!   RORR : random overlap with resorting and rebinning (Amundsen et al. 2017)
!   AEE  : adaptive equivalent extinction - the major absorber's k-distribution plus the others'
!          arithmetic band means (Amundsen et al. 2017)
!   EE   : as AEE, but each minor gas's grey value in each layer is its k-distribution weighted
!          by the flux each g-ordinate carries with that gas alone (Edwards 1996):
!          w_g ~ Gw_g exp(-tau_g/mu_z) for the stellar beam (shortwave set), and for the longwave
!          set Gw_g exp(-D tau_above) (EE_lw = 'escape') or Gw_g [exp(-D tau_above) +
!          exp(-D tau_below)] (EE_lw = 'twoside', default), D = 1.66, tau at the layer middle
!          for that gas alone (including the column above the model top). The major gas, the one
!          with the largest grey optical depth from the top to where the total reaches 1, keeps
!          its k-distribution. EE returns two opacity sets: shortwave (k_sw, a_sw, g_sw, needs
!          mu_z) and longwave (k_tot, a_tot, g_tot); the other methods return the same set twice.
!!!

module ck_opacity_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = REAL64

  !! Common constants
  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp*pi, fourpi = 4.0_dp*pi
  real(dp), parameter :: half = 1.0_dp/2.0_dp
  real(dp), parameter :: third = 1.0_dp/3.0_dp, twothird = 2.0_dp/3.0_dp

  real(dp), parameter :: amu = 1.66053906660e-24_dp
  real(dp), parameter :: amu_cgs = amu
  real(dp), parameter :: sb  = 5.670374419e-8_dp
  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: kb_cgs = kb
  real(dp), parameter :: R_gas = 8.31446261815324_dp
  real(dp), parameter :: hp = 6.62607015e-34_dp
  real(dp), parameter :: c_s = 2.99792458e8_dp
  real(dp), parameter :: kb_si = 1.380649e-23_dp

  type ck_table

    character(len=10) :: sp  ! Species name
    character(len=250) :: path ! Path to table data
    integer :: nwl
    real(dp), allocatable, dimension(:) :: wl ! Central wavelengths of table (um)
    real(dp), allocatable, dimension(:) :: wn ! Central wavenumbers of table (cm-1)

    integer :: iVMR ! VMR index for consituent species

    integer :: nT ! Number of temperature points in table
    real(dp), allocatable, dimension(:) :: T  ! Temperature points in table (K)
    real(dp), allocatable, dimension(:) :: lT  ! Temperature points in table (K)

    integer :: nP ! Number of pressure points in table
    real(dp), allocatable, dimension(:) :: P ! Pressure points in table (bar) [CARE UNITS sometimes!!]
    real(dp), allocatable, dimension(:) :: lP ! Pressure points in table (bar) [CARE UNITS sometimes!!]

    integer :: nG ! Number of g-ordinances (for corr-k)
    real(dp), allocatable, dimension(:) :: Gx, Gw ! Gx ordinances in table, Gw weight of ordinances

    real(dp), allocatable, dimension(:,:,:,:) :: kap ! Kappa values [cm2 molecule-1] (CARE: convert units from table source)

  end type ck_table

  type CIA_table

    character(len=10) :: sp  ! Species name
    logical :: i3 = .False. ! Flag for 3 part special species
    character(len=10), dimension(2) :: sp_con ! Constituent species (for CIA)
    character(len=10), dimension(3) :: sp_con_3 ! Constituent species for 3 part special
    integer, dimension(2) :: iVMR ! VMR index for consituent species
    integer, dimension(3) :: iVMR_3 ! VMR index for 3 part special
    character(len=250) :: path ! Path to table data

    integer :: form

    integer :: nwl
    real(dp), allocatable, dimension(:) :: wl ! Central wavelengths of table (um)
    real(dp), allocatable, dimension(:) :: wn ! Central wavenumbers of table (cm-1)

    integer :: nT ! Number of temperature points in table
    real(dp), allocatable, dimension(:) :: T  ! Temperature points in table (K)

    real(dp), allocatable, dimension(:,:) :: kap ! table values [Usually: cm5 molecule-2] (CARE: Check format units)

    !! Interpolation axes and log10 table values (set after reading): lx - wn (form 1) or
    !! log10 wl [A] (form 2); lT - log10 T (forms 1, 3) or log10 theta (form 2); lkap - log10 kap
    real(dp), allocatable, dimension(:) :: lx, lT
    real(dp), allocatable, dimension(:,:) :: lkap

  end type CIA_table

  type Ray_table

    character(len=10) :: sp  ! Species name
    integer :: iVMR ! VMR index for consituent species
    character(len=250) :: path ! Path to table data
    integer :: nwl
    real(dp), allocatable, dimension(:) :: kap

  end type Ray_table

  integer :: u_nml

  logical :: first_call = .True.
  integer :: n_wl, n_b, n_g
  real(dp), allocatable, dimension(:) :: wl_e, wl, wn_e, wn, freq
  character(len=250) :: wl_path

  ! CK table data - read in pre-calculated cm2 molecule-1 k-coefficents for each band
  logical :: PM, RORR, AEE, EE
  character(len=10) :: EE_lw
  integer :: ee_lw_kind
  real(dp), parameter :: D_lw = 1.66_dp               ! diffusivity factor of the EE longwave weights
  integer :: ck_interp = 2              ! 1 = (bi)linear, 2 = Bezier, for the k-tables and CIA
  integer :: n_ck
  integer, allocatable, dimension(:) :: ck_form
  type(ck_table), allocatable, dimension(:) :: ck
  character(len=10), allocatable, dimension(:) :: ck_sp
  character(len=250), allocatable, dimension(:) :: ck_paths

    ! CIA table data - read in HITRAN based CIA table then interpolate
  integer :: n_CIA
  integer, allocatable, dimension(:) :: CIA_form
  type(CIA_table), allocatable, dimension(:) :: CIA
  character(len=10), allocatable, dimension(:) :: CIA_sp
  character(len=250), allocatable, dimension(:) :: CIA_paths

  ! Rayleigh table data - read in pre-calcuilated cm2 molecule-1 xsections for each band
  integer :: n_Ray
  integer, allocatable, dimension(:) :: Ray_form
  type(Ray_table), allocatable, dimension(:) :: Ray
  character(len=10), allocatable, dimension(:) :: Ray_sp
  character(len=250), allocatable, dimension(:) :: Ray_paths

  namelist /ck_nml/ PM, RORR, AEE, EE, EE_lw, ck_sp, ck_paths, ck_form, ck_interp
  namelist /CIA_nml/ CIA_sp, CIA_paths, CIA_form
  namelist /Ray_nml/ Ray_sp, Ray_paths, Ray_form


contains

  subroutine ck_opacity(nlay, n_ck, n_CIA, n_Ray, n_b, n_g, wl_in, grav_in, T_in, p_in, pe_in, mu_in, n_sp, sp_list, VMR, &
    & k_tot, a_tot, g_tot, mu_z, k_sw, a_sw, g_sw)
    implicit none

    ! What's required as input is the layer T, P, VMR
    integer, intent(in) :: nlay, n_ck, n_CIA, n_Ray, n_sp, n_b, n_g
    character(len=10), dimension(n_sp), intent(in) :: sp_list
    real(dp), intent(in) :: grav_in
    real(dp), dimension(nlay), intent(in) :: T_in, p_in, mu_in
    real(dp), dimension(nlay+1), intent(in) :: pe_in
    real(dp), dimension(n_b+1), intent(in) :: wl_in
    real(dp), dimension(n_sp, nlay), intent(in) :: VMR
    real(dp), dimension(n_g,n_b,nlay), intent(out) :: k_tot, a_tot, g_tot

    !! Shortwave set (EE: its own minor-gas weighting with mu_z the cos(zenith angle);
    !! the other methods: equal to k_tot, a_tot, g_tot)
    real(dp), intent(in), optional :: mu_z
    real(dp), dimension(n_g,n_b,nlay), intent(out), optional :: k_sw, a_sw, g_sw

    integer :: b, g, k, s
    real(dp), allocatable, dimension(:,:,:), save :: k_ck_sp
    real(dp), allocatable, dimension(:,:), save :: k_ck
    real(dp), dimension(n_ck,n_g,nlay) :: km
    real(dp), dimension(n_g,nlay) :: k_ck_sw
    real(dp), dimension(nlay) :: k_cont, k_Ray, dm

    real(dp) :: grav, mu_sw
    logical :: want_sw
    real(dp), dimension(nlay) :: Tl, pl, Nl, RH, mu
    real(dp), dimension(nlay+1) :: pe

    logical, save :: first_call = .True.

    if (first_call .eqv. .True.) then
      call read_gas_data(n_ck, n_CIA, n_Ray, n_sp, n_b, sp_list)
      !! Read input variables from namelist
      allocate(k_ck_sp(n_ck,n_g,nlay))
      allocate(k_ck(n_g,nlay))

      ! Find wavelength centers and wavenumbers, freq etc
      allocate(wl_e(n_b+1), wl(n_b), wn_e(n_b+1), wn(n_b), freq(n_b))
      do b = 1, n_b
        wl_e(b) = wl_in(b)
        wn_e(b) = 1.0_dp/(wl_e(b) * 1.0e-4_dp)
        wl(b) = (wl_in(b+1)+wl_in(b))/2.0_dp
        wn(b) = 1.0_dp/(wl(b) * 1.0e-4_dp)
        freq(b) = c_s / (wl(b) * 1.0e-6_dp)
      end do
      wl_e(n_b+1) = wl_in(n_b+1)
      wn_e(n_b+1) = 1.0_dp/(wl_in(n_b+1) * 1.0e-4_dp)

      first_call = .False.
    end if

    ! Keep everything in CGS units - p in bar
    grav = grav_in*100.0_dp ! gravity in cm s-2
    Tl(:) = T_in(:)
    pe(:) = pe_in(:)*10.0_dp ! pe in dyne
    pl(:) = p_in(:)/1e5_dp
    mu(:) = mu_in(:)
    Nl(:) = (pl(:)*1e6_dp) / (kb_cgs * Tl(:))
    RH(:) = (pl(:)*1e6_dp * mu(:) * amu_cgs) / (kb_cgs * Tl(:))
    dm(:) = (pe(2:nlay+1) - pe(1:nlay))/grav    ! layer mass [g cm-2]

    want_sw = present(k_sw)
    if (want_sw .and. EE .and. .not. present(mu_z)) stop 'ck_opacity: EE needs mu_z for k_sw'
    mu_sw = 1.0_dp
    if (present(mu_z)) mu_sw = max(mu_z, 1.0e-3_dp)   ! night: the shortwave set is not used

    ! Loop through each band - we follow various recipes from a mix of CMCRT, CHIMERA and NASA AMES
    do b = 1, n_b

      ! First, bi-linearly interp each species k-tables to mid point T and P
      ! Then perform k-coefficent mixing (Random overlap)
      if (n_ck > 0) then

        !print*, b, 'ck_interp'
        do k = 1, nlay
          if (ck_interp == 1) then
            call interp_ck_sp_TP(n_ck, n_g, Tl(k), pl(k), b, k_ck_sp(:,:,k))
          else if (ck_interp == 2) then
            call interp_ck_sp_TP_Bezier(n_ck, n_g, Tl(k), pl(k), b, k_ck_sp(:,:,k))
          end if
        end do

        if (RORR .eqv. .True.) then
          !print*, b, 'RO'
          do k = 1, nlay
            call random_overlap(n_ck, n_g, n_sp, Nl(k), RH(k), VMR(:,k), k_ck_sp(:,:,k), k_ck(:,k))
          end do
        else if (AEE .eqv. .True.) then
          !print*, b, 'AEE'
          call adap_equiv_extinction(nlay, n_ck, n_g, n_sp, grav, &
          &  pe(:), Nl(:), RH(:), VMR(:,:), k_ck_sp(:,:,:), k_ck(:,:))
        else if (EE .eqv. .True.) then
          ! Mass opacity of each species [cm2 g-1]
          do k = 1, nlay
            do s = 1, n_ck
              km(s,:,k) = k_ck_sp(s,:,k) * VMR(ck(s)%iVMR,k) * Nl(k) / RH(k)
            end do
          end do
          call flux_equiv_extinction(nlay, n_ck, n_g, km, dm, pe(1)/grav, ee_lw_kind, 0.0_dp, k_ck)
          if (want_sw) call flux_equiv_extinction(nlay, n_ck, n_g, km, dm, pe(1)/grav, 0, mu_sw, k_ck_sw)
        else if (PM .eqv. .True.) then
          !print*, b, 'PM'
          k_ck(:,:) = k_ck_sp(1,:,:)
        end if
      else
        k_ck(:,:) = 0.0_dp
      end if
      if (.not. ((EE .eqv. .True.) .and. (n_ck > 0))) k_ck_sw(:,:) = k_ck(:,:)

      ! Second, interp to find  continuum opacity
      if (n_CIA > 0) then
        !print*, b, 'conti'
        do k = 1, nlay
          if (ck_interp == 1) then
            call interp_conti(n_CIA, n_sp, n_b, b, Tl(k), Nl(k), RH(k), VMR(:,k), k_cont(k))
          else
            call interp_conti_Bezier(n_CIA, n_sp, n_b, b, Tl(k), Nl(k), RH(k), VMR(:,k), k_cont(k))
          end if
        end do
      else
        k_cont(:) = 0.0_dp
      end if

      ! Third, calculate Rayleigh scattering
      if (n_Ray > 0) then
        !print*, b, 'Ray'
        do k = 1, nlay
          call calc_Rayleigh(n_Ray, n_sp, b, Nl(k), RH(k), VMR(:,k), k_Ray(k))
        end do
      else
        k_Ray(:) = 0.0_dp
      end if

      ! Fourth, aerosol opacity
      ! (passed to subroutine directly)

      ! Find total values for each g ordinate
      do g = 1, n_g
         k_tot(g,b,:) = k_ck(g,:) + k_cont(:) + k_Ray(:) !+ k_cl_ext(b)  ! Total extinction
         a_tot(g,b,:) = min(k_Ray(:)/k_tot(g,b,:),0.99_dp) !(k_Ray + k_cl_sca(b)) / k_tot(b,g)           ! Effective single scattering albedo
         g_tot(g,b,:) = 0.0_dp !(g_cl(b)*k_cl_sca(b)) / (k_cl_sca(b) + k_Ray) ! Effective asymmetry factor
         if (want_sw) then
           k_sw(g,b,:) = k_ck_sw(g,:) + k_cont(:) + k_Ray(:)
           a_sw(g,b,:) = min(k_Ray(:)/k_sw(g,b,:),0.99_dp)
           g_sw(g,b,:) = g_tot(g,b,:)
         end if
      end do
      !print*, b, k_ck(1), k_ck(8), k_Ray, k_cont, a_tot(b,1) , a_tot(b,8)
    end do

    ! MKS conversion
    k_tot(:,:,:) = k_tot(:,:,:) * 0.1_dp
    if (want_sw) k_sw(:,:,:) = k_sw(:,:,:) * 0.1_dp

    !stop

  end subroutine ck_opacity

  !! Flux-weighted equivalent extinction of one band and column (EE). km(n_ck, n_g, nlay)
  !! [cm2/g] each gas's mass opacity, dm(nlay) layer mass [g/cm2], m_top the mass above the
  !! top level. wkind 0: stellar beam at mu (w ~ Gw exp(-tau/mu)); 1: longwave escape
  !! (Gw exp(-D tau_above)); 2: longwave two-sided (Gw [exp(-D tau_above) + exp(-D tau_below)]);
  !! tau at the layer middle for that gas alone. Returns k_ck(n_g, nlay) [cm2/g]: the major
  !! gas's k-distribution plus the others' weighted grey values.
  subroutine flux_equiv_extinction(nlay, n_ck, n_g, km, dm, m_top, wkind, mu, k_ck)
    implicit none

    integer, intent(in) :: nlay, n_ck, n_g, wkind
    real(dp), dimension(n_ck,n_g,nlay), intent(in) :: km
    real(dp), dimension(nlay), intent(in) :: dm
    real(dp), intent(in) :: m_top, mu
    real(dp), dimension(n_g,nlay), intent(out) :: k_ck

    integer :: s, k, imax
    real(dp), dimension(n_g) :: tau_run, tau_col, tmid, lw, w, a, c
    real(dp), dimension(n_ck,nlay) :: ke
    real(dp), dimension(n_ck) :: tau_e
    real(dp) :: tau_tot

    do s = 1, n_ck
      tau_col(:) = km(s,:,1)*m_top
      do k = 1, nlay
        tau_col(:) = tau_col(:) + km(s,:,k)*dm(k)
      end do
      tau_run(:) = km(s,:,1)*m_top
      do k = 1, nlay
        tmid(:) = tau_run(:) + 0.5_dp*km(s,:,k)*dm(k)
        select case (wkind)
        case (0)
          lw(:) = -tmid(:)/mu
        case (1)
          lw(:) = -D_lw*tmid(:)
        case default
          a(:) = -D_lw*tmid(:)
          c(:) = -D_lw*(tau_col(:) - tmid(:))
          lw(:) = max(a(:), c(:)) + log(1.0_dp + exp(min(a(:), c(:)) - max(a(:), c(:))))
        end select
        w(:) = ck(s)%Gw(:) * exp(lw(:) - maxval(lw(:)))
        ke(s,k) = sum(w(:)*km(s,:,k)) / sum(w(:))
        tau_run(:) = tau_run(:) + km(s,:,k)*dm(k)
      end do
    end do

    !! The major gas: largest grey optical depth from the top to where the total reaches 1
    tau_e(:) = ke(:,1)*m_top
    tau_tot = sum(tau_e(:))
    do k = 1, nlay
      if (tau_tot >= 1.0_dp) exit
      tau_e(:) = tau_e(:) + ke(:,k)*dm(k)
      tau_tot = sum(tau_e(:))
    end do
    imax = maxloc(tau_e, dim=1)

    do k = 1, nlay
      k_ck(:,k) = km(imax,:,k) + (sum(ke(:,k)) - ke(imax,k))
    end do

  end subroutine flux_equiv_extinction

  subroutine adap_equiv_extinction(nlay, n_ck, n_g, n_sp, grav, pe, Nl, RH, VMR, k_ck_sp, k_ck)
    implicit none

    integer, intent(in) :: nlay, n_ck, n_g, n_sp
    real(dp), intent(in) :: grav
    real(dp), dimension(nlay), intent(in) :: Nl, RH
    real(dp), dimension(nlay+1), intent(in) :: pe
    real(dp), dimension(n_sp,nlay), intent(in) :: VMR
    real(dp), dimension(n_ck,n_g,nlay), intent(in) :: k_ck_sp

    real(dp), dimension(n_g,nlay), intent(out) :: k_ck

    integer :: k, s, g
    integer, dimension(n_ck) :: max_idx
    real(dp) :: top, bot, dpe, tau_tot, tau_s, tau_all
    real(dp), dimension(n_ck,nlay) :: k_av
    real(dp), dimension(n_ck) :: tau_av
    logical, dimension(n_ck) :: major
    !real(dp), dimension(n_g) :: Trans, tau_s_g

    !! Calculate average kappa in for each species in each layer for thermal component
    !! First just try a general non-weighted average
    !! - seems a good compromise for both stellar and thermal components
    do s = 1, n_ck
      do k = 1, nlay
        !bl(k) = BB(wl_in, Tl(k))
        top = 0.0_dp
        bot = 0.0_dp
        do g = 1, n_g
          top = top + ck(s)%Gw(g)*k_ck_sp(s,g,k) !*bl(k)
          bot = bot + ck(s)%Gw(g) !*bl(k)
        end do
        k_av(s,k) = top/bot * (VMR(ck(s)%iVMR,k) * Nl(k) / RH(k))
      end do
    end do

    !! Find the major absorber species
    !! If we use mean opacity for all species, then using tau is fine
    !! no need to do the transmission functions
    tau_av(:) = 0.0_dp
    tau_tot = 0.0_dp
    do k = 1, nlay
      dpe = pe(k+1) - pe(k)
      do s = 1, n_ck
        tau_s = (k_av(s,k) * dpe)/grav
        tau_av(s) = tau_av(s) + tau_s
      end do
      tau_all = (sum(k_av(:,k)) * dpe)/grav
      tau_tot = tau_tot + tau_all
      !print*, k, pe(k), tau_av(:), tau_tot
      if (tau_tot >= 1.0_dp) then
        exit
      end if
    end do

    !! Find species with highest tau_av
    !! This also checks for the largest tau should the surface be at tau_av < 1
    max_idx = maxloc(tau_av,dim=1)
    !print*, max_idx
    major(:) = .False.
    major(max_idx(1)) = .True.

    !! Combine the grey average opacities of minor species with the k-coefficents of the major species
    k_ck(:,:) = 0.0_dp
    do k = 1, nlay
     do s = 1, n_ck
       if (major(s) .eqv. .True.) then
         !! Combine full k-table
         k_ck(:,k) = k_ck(:,k) + k_ck_sp(s,:,k) * (VMR(ck(s)%iVMR,k) * Nl(k) / RH(k))
       else
         !! Combine average grey value
         k_ck(:,k) = k_ck(:,k) + k_av(s,k)
       end if
     end do
    end do

    !! Combine the grey average opacities of minor species with the k-coefficents of the major species
    !! Multiply the transmission of each species at each layer for each g coordinate
    ! do k = 1, nlay
    !   Trans(:) = 1.0_dp
    !   dpe = pe(k+1) - pe(k)
    !   do s = 1, n_ck
    !     if (major(s) .eqv. .True.) then
    !       ! Combine full k-table
    !       tau_s_g(:) = ((k_ck_sp(s,:,k) * VMR(ck(s)%iVMR,k) * Nl(k) / RH(k)) *  dpe)/grav
    !       Trans(:) =  Trans(:) * exp(-tau_s_g(:))
    !     else
    !       ! Combine average grey value
    !       tau_s = (k_av(s,k) * dpe)/grav
    !       Trans(:) = Trans(:) * exp(-tau_s)
    !     end if
    !   end do
    !   ! Find the opacity of each g-ordinance for this layer
    !   k_ck(:,k) = -(log(Trans(:)) * grav) / dpe
    !   print*,k_ck(:,k) 
    ! end do


  end subroutine adap_equiv_extinction

  subroutine random_overlap(n_ck, n_g,  n_sp, N, RH, VMR, k_ck_sp, k_ck)
    implicit none

    integer, intent(in) :: n_ck, n_g, n_sp
    real(dp), intent(in) :: N, RH
    real(dp), dimension(n_sp), intent(in) :: VMR
    real(dp), dimension(n_ck,n_g), intent(in) :: k_ck_sp
    real(dp), dimension(n_g), intent(out) :: k_ck

    integer :: i, j, g, s, n_g2
    real(dp) :: VMR_tot, VMR_cum
    real(dp), dimension(n_g*n_g) :: k_mix
    real(dp), dimension(n_g*n_g) :: wt_mix
    real(dp), dimension(0:n_g*n_g) :: intg, x

    integer :: ix, ix1
    real(dp) :: xval

    !! If 1 ck table, then it's just VMR * ck coefficents
    if (n_ck == 1)then
      k_ck(:) =  VMR(ck(1)%iVMR) * k_ck_sp(1,:) * N / RH
      return
    end if

    !! Convolved size
    n_g2 = n_g * n_g

    !! Start RO procedure
    ! Initial k-table is 1st table * VMR
    VMR_tot = VMR(ck(1)%iVMR)
    k_ck(:) = k_ck_sp(1,:)
    ! Loop over all other tables
    do s = 2, n_ck
      ! Track current cumulative VMR (_cum) and next VMR (_tot)
      VMR_cum = VMR_tot
      VMR_tot = VMR_tot + VMR(ck(s)%iVMR)
      ! Skip this species if very low abundance
      if (VMR(ck(s)%iVMR)*sum(k_ck_sp(s,:)*ck(s)%Gw(:)) < 1.0e-99_dp) then
        cycle
      end if
      ! Loop n_g by n_g, perform random k and weight mixing following Amundsen et al. (2017)
      do i = 1, n_g
        do j = 1, n_g
          k_mix((i-1)*n_g+j) = (VMR_cum*k_ck(i) + VMR(ck(s)%iVMR)*k_ck_sp(s,j)) &
            & / VMR_tot
          wt_mix((i-1)*n_g+j) = ck(s-1)%Gw(i) * ck(s)%Gw(j)
        end do
      end do

      ! Sort the mixed k tables and assosiated weights
      call sort2(n_g2, k_mix, wt_mix)

      ! Now reconstruct the x (g) coordinate
      ! Find cumulative sum of the mixed weights
      intg(0) = 0.0_dp
      do g = 1, n_g2
        intg(g) = intg(g-1) + wt_mix(g)
      end do
      ! Each mixed k at the centre of its normalised cumulative-weight bin
      x(1:n_g2) = (intg(0:n_g2-1) + 0.5_dp*wt_mix(:))/intg(n_g2)

      ! Now interpolate to the origional x grid, this is the mixed k-table.
      ! (Exo-FMS_2 fix: the old search ran through x(0:n) via a view starting at 1, so it
      ! interpolated between the wrong pair of cumulative weights, one bin too high)
      do g = 1, n_g
        ! x(ix) < Gx <= x(ix+1); the end values outside
        xval = ck(s)%Gx(g)
        call locate(x(1:n_g2), xval, ix)
        if (ix < 1) then
          k_ck(g) = k_mix(1)
        else if (ix >= n_g2) then
          k_ck(g) = k_mix(n_g2)
        else
          ix1 = ix + 1
          call linear_log_interp(xval, x(ix), x(ix1), k_mix(ix), k_mix(ix1), k_ck(g))
        end if
      end do

    end do

    k_ck(:) = VMR_tot * N * k_ck(:) / RH

  end subroutine random_overlap

  subroutine interp_ck_sp_TP(n_ck,n_g, T, P, b, k_ck_sp)
    implicit none

    integer, intent(in) :: n_ck, n_g, b
    real(dp), intent(in) :: T, P
    real(dp), dimension(n_ck,n_g), intent(out) :: k_ck_sp

    integer :: s, g
    integer :: iT, iT1, iP, iP1
    real(dp) :: lT, lP
    real(dp) :: xval, x0, x1, y0, y1, yval, aval
    real(dp) :: a00, a01, a10, a11

    !! A large if block statement is used to determine if the T-p point
    !! is within the ck table, then interpolates.
    !! If T > Tmax, then T = Tmax, if P > Pmax, P = Pmax and if P < Pmin, P = Pmin
    !! This is done to avoid crashes with k-tables that may not cover the intendend T-p range
    !! IT IS AVDISED TO CREATE k-tables WITHIN APPRORIATE T-p RANGES!!
    !! NOTE: no such check is done for T < Tmin of the table, which will give an error

    lT = log10(T)
    lP = log10(P)

    do s = 1, n_ck

      ! Find temperature grid index
      call locate(ck(s)%T(:),T,iT)
      iT1 = iT + 1
      ! Find pressure grid index
      call locate(ck(s)%P(:),P,iP)
      iP1 = iP + 1

      if (iT == ck(s)%nT) then
        ! Temperature is too high, outside table range - use highest availible T data
        if (iP == ck(s)%nP) then
           ! Pressure is too high, outside table range - use highest availible P data
            k_ck_sp(s,:) = 10.0_dp**ck(s)%kap(b,iT,iP,:)
        else if (iP == 0) then
           ! Pressure is too low, outside table range - use highest lowest P data
           k_ck_sp(s,:) = 10.0_dp**ck(s)%kap(b,iT,1,:)
        else
          ! Pressure is within table range, perform linear interpolation at highest T
          do g = 1, ck(s)%nG
            xval = lP
            x0 = ck(s)%lP(iP) ; x1 = ck(s)%lP(iP1)
            y0 = ck(s)%kap(b,iT,iP,g) ; y1 = ck(s)%kap(b,iT,iP1,g)
            call linear_interp(xval, x0, x1, y0, y1, yval)
            !call linear_log_interp(xval, x0, x1, y0, y1, yval)
            k_ck_sp(s,g) = 10.0_dp**yval
          end do
        end if
      else if (iT == 0) then
        ! Temperature is too low, outside table range - use lowest availible T data
        if (iP == ck(s)%nP) then
           ! Pressure is too high, outside table range - use highest availible P data
            k_ck_sp(s,:) = 10.0_dp**ck(s)%kap(b,1,iP,:)
        else if (iP == 0) then
           ! Pressure is too low, outside table range - use highest lowest P data
           k_ck_sp(s,:) = 10.0_dp**ck(s)%kap(b,1,1,:)
        else
          ! Pressure is within table range, perform linear interpolation at highest T
          do g = 1, ck(s)%nG
            xval = lP
            x0 = ck(s)%lP(iP) ; x1 = ck(s)%lP(iP1)
            y0 = ck(s)%kap(b,1,iP,g) ; y1 = ck(s)%kap(b,1,iP1,g)
            call linear_interp(xval, x0, x1, y0, y1, yval)
            !call linear_log_interp(xval, x0, x1, y0, y1, yval)
            k_ck_sp(s,g) = 10.0_dp**yval
          end do
        end if
      else
        ! Temperature is within the normal range
        if (iP == ck(s)%nP) then
          ! Pressure is too high, outside table range - use highest availible P
          do g = 1, ck(s)%nG
            xval = lT
            x0 = ck(s)%lT(iT) ; x1 = ck(s)%lT(iT1)
            y0 = ck(s)%kap(b,iT,iP,g) ; y1 = ck(s)%kap(b,iT1,iP,g)
            call linear_interp(xval, x0, x1, y0, y1, yval)
            !call linear_log_interp(xval, x0, x1, y0, y1, yval)
            k_ck_sp(s,g) = 10.0_dp**yval
          end do
        else if (iP == 0) then
          ! Pressure is too low, outside table range - use lowest availible P
          do g = 1, ck(s)%nG
            xval = lT
            x0 = ck(s)%lT(iT) ; x1 = ck(s)%lT(iT1)
            y0 = ck(s)%kap(b,iT,1,g) ; y1 = ck(s)%kap(b,iT1,1,g)
            call linear_interp(xval, x0, x1, y0, y1, yval)
            !call linear_log_interp(xval, x0, x1, y0, y1, yval)
            k_ck_sp(s,g) = 10.0_dp**yval
          end do
        else
          ! Both pressure and temperature are within table bounds, perform bi-linear interpolation
          xval = lT
          yval = lP
          x0 = ck(s)%lT(iT) ; x1 = ck(s)%lT(iT1)
          y0 = ck(s)%lP(iP) ; y1 = ck(s)%lP(iP1)
          do g = 1, ck(s)%nG
            a00 = ck(s)%kap(b,iT,iP,g) ; a01 = ck(s)%kap(b,iT,iP1,g)
            a10 = ck(s)%kap(b,iT1,iP,g) ; a11 = ck(s)%kap(b,iT1,iP1,g)
            call bilinear_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)
            !call bilinear_log_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)
            k_ck_sp(s,g) = 10.0_dp**aval
          end do
        end if
      end if

    end do

  end subroutine interp_ck_sp_TP

  subroutine interp_ck_sp_TP_Bezier(n_ck,n_g, T, P, b, k_ck_sp)
    implicit none

    integer, intent(in) :: n_ck, n_g, b
    real(dp), intent(in) :: T, P
    real(dp), dimension(n_ck,n_g), intent(out) :: k_ck_sp

    integer :: s, g, iT, iP
    logical :: eT, eP
    real(dp) :: lT, lP

    !! Bezier interpolation of log10 k in log10 T and log10 p (shared with Exo-FMS_2), table
    !! edge values outside the table.
    !! IT IS AVDISED TO CREATE k-tables WITHIN APPRORIATE T-p RANGES!!

    lT = log10(T)
    lP = log10(P)

    do s = 1, n_ck
      call bezier_index(ck(s)%lT(:), lT, iT, eT)
      call bezier_index(ck(s)%lP(:), lP, iP, eP)
      do g = 1, n_g
        k_ck_sp(s,g) = 10.0_dp**bezier_2D(ck(s)%lT(:), ck(s)%lP(:), ck(s)%kap(b,:,:,g), lT, lP, iT, eT, iP, eP)
      end do
    end do

  end subroutine interp_ck_sp_TP_Bezier

  subroutine interp_conti(n_CIA, n_sp, n_b, b,  T, N, RH, VMR, k_cont)
    implicit none

    integer, intent(in) :: b, n_CIA, n_sp, n_b
    real(dp), intent(in) :: T, N, RH
    real(dp), dimension(n_sp), intent(in) :: VMR
    real(dp), intent(out) :: k_cont

    integer, allocatable, dimension(:,:), save :: iwn, iwn1
    logical, save :: f_call = .True.

    integer :: s
    integer :: iT, iT1
    real(dp) :: k_Hm, k_ff
    real(dp) :: x0, x1, y0, y1, xval, yval
    real(dp) :: a00, a01, a10, a11, aval

    ! If this is the first call, find the wavenumber index numbers
    ! so we don't have to calculate them each time
    if (f_call .eqv. .True.) then
      if (allocated(iwn) .eqv. .False.) then
        allocate(iwn(n_CIA,n_b), iwn1(n_CIA,n_b))
      end if
      do s = 1, n_CIA
        if ((trim(CIA(s)%sp) == 'H-') .or. (CIA(s)%form == 2)) then
          cycle
        end if
        call locate(CIA(s)%wn(:),wn(b),iwn(s,b))
        iwn1(s,b) = iwn(s,b) + 1
        !print*, wn(b), iwn(s,b), iwn1(s,b), CIA(s)%wn(iwn(s,b)), CIA(s)%wn(iwn1(s,b))
      end do
      if (b == n_b) then
        f_call = .False.
      end if
    end if

    k_cont = 0.0_dp
    do s = 1, n_CIA

      if (trim(CIA(s)%sp) == 'H-') then
        call CIA_Hminus(n_sp, s, b, T, N, VMR, k_Hm)
        k_cont = k_cont + k_Hm
        cycle
      end if
      if (CIA(s)%form == 2) then
        call CIA_ff(n_sp, s, b, T, N, VMR, k_ff)
        k_cont = k_cont + k_ff
        !print*, s, b, k_ff
        cycle
      end if
      if (CIA(s)%form == 3) then
        ! Table already on the model's bands - linear interpolation in T (as Exo-FMS_2),
        ! table edge values outside the T range
        k_cont = k_cont + CIA_T_linear(s, b, T) &
          & * VMR(CIA(s)%iVMR(1)) * N &
          & * VMR(CIA(s)%iVMR(2)) * N
        cycle
      end if

      ! If band is outside table wavelength range
      if ((iwn1(s,b) > CIA(s)%nwl) .or. (iwn(s,b) < 1)) then
        cycle
      end if

      ! Locate required T indexes in CIA wn array for layer temperature
      call locate(CIA(s)%T(:),T,iT)
      iT1 = iT + 1

      !! Perform temperature edge case check
      if (iT < 1) then
        ! Temperature of layer is outside lower bounds of table
        ! Perform wn linear interp to minval(T)
        xval = wn(b) ; x0 = CIA(s)%wn(iwn(s,b)) ; x1 = CIA(s)%wn(iwn1(s,b))
        y0 = CIA(s)%kap(iwn(s,b),1) ; y1 = CIA(s)%kap(iwn1(s,b),1)

        ! Perform log linear interpolation
        call linear_interp(xval, x0, x1, y0, y1, yval)

        ! Add to result to work variable in units of [cm-1]
        k_cont = k_cont + yval &
          & * VMR(CIA(s)%iVMR(1)) * N &
          & * VMR(CIA(s)%iVMR(2)) * N

      else if (iT1 > CIA(s)%nT) then

        ! Temperature of layer is outside upper bounds of table
        ! Perform wn linear interp to maxval(T)
        xval = wn(b) ; x0 = CIA(s)%wn(iwn(s,b)) ; x1 = CIA(s)%wn(iwn1(s,b))
        y0 = CIA(s)%kap(iwn(s,b),CIA(s)%nT) ; y1 = CIA(s)%kap(iwn1(s,b),CIA(s)%nT)

        ! Perform log linear interpolation
        call linear_interp(xval, x0, x1, y0, y1, yval)

        ! Add to result to work variable in units of [cm-1]
        k_cont = k_cont + yval &
          & * VMR(CIA(s)%iVMR(1)) * N &
          & * VMR(CIA(s)%iVMR(2)) * N

      else

        !! wn and T are within the table bounds
        xval = wn(b) ; x0 = CIA(s)%wn(iwn(s,b)) ; x1 = CIA(s)%wn(iwn1(s,b))
        yval = T ; y0 = CIA(s)%T(iT) ; y1 = CIA(s)%T(iT1)
        a00 = CIA(s)%kap(iwn(s,b),iT) ; a10 = CIA(s)%kap(iwn1(s,b),iT)
        a01 = CIA(s)%kap(iwn(s,b),iT1) ; a11 = CIA(s)%kap(iwn1(s,b),iT1)

        ! Perform bi-linear interpolation
        call bilinear_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)

        ! Add to result to work variable in units of [cm-1]
        k_cont = k_cont + aval &
          & * VMR(CIA(s)%iVMR(1)) * N &
          & * VMR(CIA(s)%iVMR(2)) * N

      end if

    end do

    k_cont = k_cont / RH

  end subroutine interp_conti


  subroutine interp_conti_Bezier(n_CIA, n_sp, n_b, b,  T, N, RH, VMR, k_cont)
    implicit none

    integer, intent(in) :: b, n_CIA, n_sp, n_b
    real(dp), intent(in) :: T, N, RH
    real(dp), dimension(n_sp), intent(in) :: VMR
    real(dp), intent(out) :: k_cont

    integer :: s, ix, iT
    logical :: ex, eT
    real(dp) :: k_Hm, k_ff, kap, lT

    !! Bezier interpolation (shared with Exo-FMS_2) of log10 kap, table edge values outside
    !! the table (form 1: zero outside the table's wavenumber range, where there is no data)

    lT = log10(T)

    k_cont = 0.0_dp
    do s = 1, n_CIA

      if (trim(CIA(s)%sp) == 'H-') then
        call CIA_Hminus(n_sp, s, b, T, N, VMR, k_Hm)
        k_cont = k_cont + k_Hm
        cycle
      end if

      select case (CIA(s)%form)
      case (2)
        ! Free-free table in log10 wl and log10 theta
        call CIA_ff_Bezier(n_sp, s, b, T, N, VMR, k_ff)
        k_cont = k_cont + k_ff
        cycle
      case (3)
        ! Table already on the model's bands - log10 T only
        kap = CIA_T_Bezier(s, b, T)
      case default
        ! HITRAN table in wn and log10 T, sampled at the band centre
        if ((wn(b) < CIA(s)%wn(1)) .or. (wn(b) > CIA(s)%wn(CIA(s)%nwl))) cycle
        call bezier_index(CIA(s)%lx(:), wn(b), ix, ex)
        call bezier_index(CIA(s)%lT(:), lT, iT, eT)
        kap = 10.0_dp**bezier_2D(CIA(s)%lx(:), CIA(s)%lT(:), CIA(s)%lkap(:,:), wn(b), lT, ix, ex, iT, eT)
      end select

      ! Add to result to work variable in units of [cm-1]
      k_cont = k_cont + kap &
        & * VMR(CIA(s)%iVMR(1)) * N &
        & * VMR(CIA(s)%iVMR(2)) * N

    end do

    k_cont = k_cont / RH

  end subroutine interp_conti_Bezier

  !! Form 3 CIA [cm5 molecule-2] of band b at temperature T: linear in T, edge values outside
  real(dp) function CIA_T_linear(s, b, T) result(kap)
    implicit none

    integer, intent(in) :: s, b
    real(dp), intent(in) :: T

    integer :: iT

    call locate(CIA(s)%T(:), T, iT)
    if (iT < 1) then
      kap = CIA(s)%kap(b,1)
    else if (iT + 1 > CIA(s)%nT) then
      kap = CIA(s)%kap(b,CIA(s)%nT)
    else
      call linear_interp(T, CIA(s)%T(iT), CIA(s)%T(iT+1), CIA(s)%kap(b,iT), CIA(s)%kap(b,iT+1), kap)
    end if

  end function CIA_T_linear

  !! Form 3 CIA [cm5 molecule-2] of band b at temperature T: Bezier interpolation of log10 k
  !! in log10 T (linear if the table has fewer than 3 temperatures), edge values outside
  real(dp) function CIA_T_Bezier(s, b, T) result(kap)
    implicit none

    integer, intent(in) :: s, b
    real(dp), intent(in) :: T

    kap = 10.0_dp**bezier_1D(CIA(s)%lT(:), CIA(s)%lkap(b,:), log10(T))

  end function CIA_T_Bezier

  subroutine calc_Rayleigh(n_Ray, n_sp, b, Nl, RH, VMR, k_Ray)
    implicit none

    integer, intent(in) :: b, n_Ray, n_sp
    real(dp), intent(in) :: Nl, RH
    real(dp), dimension(n_sp), intent(in) :: VMR
    real(dp), intent(out) :: k_Ray

    integer :: s

    ! For Rayleigh scattering we just need the pre-tabulated cross section
    ! in cm2 molecule-1 for each Rayleigh scattering species in this band.
    ! The only process that affects the strength of Rayleigh scattering is
    ! changes in the VMR and local total number density

    k_Ray = 0.0_dp
    do s = 1, n_Ray
      k_Ray = k_Ray + VMR(Ray(s)%iVMR) * Ray(s)%kap(b)
      !print*, s, b, k_Ray, VMR(Ray(s)%iVMR), Ray(s)%kap(b)
    end do
    k_Ray = k_Ray * Nl / RH

  end subroutine calc_Rayleigh

  ! Start up module to read in tables and prepare opacities
  subroutine read_gas_data(n_ck, n_CIA, n_Ray, n_sp, n_b,  sp_list)
    implicit none

    integer, intent(in) :: n_ck, n_CIA, n_Ray, n_sp, n_b
    character(len=10), dimension(n_sp), intent(in) :: sp_list

    integer :: u, l, i, j, g, b, n, s, uR
    integer :: ni

    real(dp), allocatable, dimension(:) :: wl_dum

    integer :: stat

    character(len=10) :: name
    integer :: nrec
    real(dp) :: wn_s, wn_e, temp_r, kmax, dum
    real(dp), allocatable, dimension(:) :: dum_b
    
    logical :: exists

    ! Open the namelist file
    open(newunit=u_nml, file='FMS_RC.nml', action='read')


    if (n_ck > 0) then

      print*, 'Reading in ck data'
      !! read in ck data and prep
      allocate(ck_sp(n_ck), ck_paths(n_ck), ck(n_ck), ck_form(n_ck))
      ! Read ck namelist variables
      PM = .False.
      RORR = .False.
      AEE = .False.
      EE = .False.
      EE_lw = 'twoside'
      ck_interp = 2
      read(u_nml, nml=ck_nml)
      if (count([PM, RORR, AEE, EE]) /= 1) then
        print*, 'ck_nml: choose exactly one of PM, RORR, AEE, EE'
        stop
      end if
      select case (trim(EE_lw))
      case ('escape')
        ee_lw_kind = 1
      case ('twoside')
        ee_lw_kind = 2
      case default
        print*, 'ck_nml: invalid EE_lw (escape or twoside): ', trim(EE_lw)
        stop
      end select
      ck(:)%sp = ck_sp(:)
      ck(:)%path = ck_paths(:)
      ! Find index for VMR
      do s = 1, n_ck
        if (PM .eqv. .True.) then
          cycle
        end if
        exists = .False.
        do j = 1, n_sp
          if (trim(ck(s)%sp) == trim(sp_list(j))) then
            ck(s)%iVMR = j
            exists = .True.
            exit
          end if
        end do
        if (exists .eqv. .False.) then
          print*, 'ERROR - Specifed CK species in namelist not found in species list - STOPPING'
          print*, 'Species: ', ck(s)%sp
          print*, 'List: ', sp_list(:)
          stop
        end if
      end do

      ! Now read in prepared ck_opacity tables
      do s = 1, n_ck
        if (ck_form(s) == 1) then
          print*, 'Reading M. Line ck table: ', trim(ck(s)%sp), ' @ ',  trim(ck(s)%path)
          open(newunit=u, file=trim(ck(s)%path), action='read')
          read(u,*) ck(s)%nwl, ck(s)%nT, ck(s)%nP, ck(s)%nG
          allocate(ck(s)%wl(ck(s)%nwl), ck(s)%wn(ck(s)%nwl))
          allocate(ck(s)%T(ck(s)%nT))
          allocate(ck(s)%P(ck(s)%nP))
          allocate(ck(s)%Gx(ck(s)%nG), ck(s)%Gw(ck(s)%nG))
          allocate(wl_dum(ck(s)%nwl+1))
          do l = 1, ck(s)%nwl+1
            read(u,*) wl_dum(l)
          end do
          do l = 1, ck(s)%nwl
             ck(s)%wl(l) = (wl_dum(l) + wl_dum(l+1)) / 2.0_dp
             ck(s)%wn(l) = 1.0_dp/ (ck(s)%wl(l) * 1e-4_dp)
             !print*, l, ck(s)%wl(l), ck(s)%wn(l)
          end do
          deallocate(wl_dum)
          do g = 1, ck(s)%nG
            read(u,*) ck(s)%Gx(g), ck(s)%Gw(g)
            !print*, g, ck(s)%Gx(g), ck(s)%Gw(g)
          end do
          do i = 1, ck(s)%nT
            read(u,*) ck(s)%T(i)
            !print*, i, ck(s)%T(i)
          end do
          do j = 1, ck(s)%nP
            read(u,*) ck(s)%P(j)
            !print*, j, ck(s)%P(j)
          end do
          read(u,*)
          allocate(ck(s)%kap(ck(s)%nwl,ck(s)%nT,ck(s)%nP,ck(s)%nG))
          do b = 1, ck(s)%nwl
            do i = 1, ck(s)%nT
              do j = 1, ck(s)%nP
                read(u,*) (ck(s)%kap(b,i,j,g), g = 1, ck(s)%nG)
                !print*, b,i,j, (ck(s)%kap(b,i,j,g), g = 1, ck(s)%nG)
              end do
            end do
          end do
        else if (ck_form(s) == 2) then
          print*, 'Reading HELIOS ck table: ', trim(ck(s)%sp), ' @ ',  trim(ck(s)%path)
          open(newunit=u, file=trim(ck(s)%path), action='read')
          read(u,*)
          read(u,*) ck(s)%nT, ck(s)%nP, ck(s)%nwl, ck(s)%nG
          allocate(ck(s)%wl(ck(s)%nwl+1), ck(s)%wn(ck(s)%nwl+1))
          allocate(ck(s)%T(ck(s)%nT),ck(s)%lT(ck(s)%nT))
          allocate(ck(s)%P(ck(s)%nP),ck(s)%lP(ck(s)%nP))
          allocate(ck(s)%Gx(ck(s)%nG), ck(s)%Gw(ck(s)%nG))

          if ((ck_interp == 2) .and. ((ck(s)%nT < 3) .or. (ck(s)%nP < 3))) then
            print*, 'Bezier k-table interpolation (ck_interp = 2) needs >= 3 T and P points: ', trim(ck(s)%sp)
            stop
          end if
          read(u,*) ck(s)%T(:)
            !print*, i, ck(s)%T(:)
          ck(s)%lT(:) = log10(ck(s)%T(:))
          read(u,*) ck(s)%P(:)
          !print*, j, ck(s)%P(:)
          ck(s)%lP(:) = log10(ck(s)%P(:))
          read(u,*)
          read(u,*) (ck(s)%wl(g), g = 1, ck(s)%nwl+1)
          read(u,*) (ck(s)%wn(g), g = 1, ck(s)%nwl+1)
          read(u,*)
          read(u,*) (ck(s)%Gx(g), g = 1, ck(s)%nG)
          read(u,*) (ck(s)%Gw(g), g = 1, ck(s)%nG)
          read(u,*)
          allocate(ck(s)%kap(ck(s)%nwl,ck(s)%nT,ck(s)%nP,ck(s)%nG))
          ! Some slight differences between premixed and not premixed opacitiy tables here
          ! Due to GGChem grid ordering in the premixed opacitiy calculation
          if (PM .eqv. .False.) then
            do i = 1, ck(s)%nT
              do j = 1, ck(s)%nP
                do b = ck(s)%nwl, 1, -1
                  read(u,*) (ck(s)%kap(b,i,j,g), g = 1, ck(s)%nG)
                  ck(s)%kap(b,i,j,:) = max(ck(s)%kap(b,i,j,:),1.0e-99_dp)
                  ck(s)%kap(b,i,j,:) = log10(ck(s)%kap(b,i,j,:))
                  !print*, b,i,j, (ck(s)%kap(b,i,j,g), g = 1, ck(s)%nG)
                end do
              end do
            end do
          else if (PM .eqv. .True.) then
            do i = 1, ck(s)%nT
              do j = 1, ck(s)%nP
                do b = ck(s)%nwl, 1, -1
                  read(u,*) (ck(s)%kap(b,i,j,g), g = 1, ck(s)%nG)
                  ck(s)%kap(b,i,j,:) = max(ck(s)%kap(b,i,j,:),1.0e-99_dp)
                  ck(s)%kap(b,i,j,:) = log10(ck(s)%kap(b,i,j,:))
                  !print*, b,i,j, (ck(s)%kap(b,j,i,g), g = 1, ck(s)%nG)
                end do
              end do
            end do
          end if
        else
          print*, 'No ck species, n_ck = ', n_ck
          stop
        end if
      end do
    end if

    if (n_CIA > 0) then

      !! read in ck data and prep
      allocate(CIA_sp(n_CIA), CIA_paths(n_CIA), CIA(n_CIA), CIA_form(n_CIA))

      ! Read the cia namelist
      read(u_nml, nml=cia_nml)

      CIA(:)%sp = CIA_sp(:)
      CIA(:)%path = CIA_paths(:)
      CIA(:)%form = CIA_form(:)

      ! Find the CIA constituents from lookup table
      call find_CIA_consituents(n_CIA)

      do s = 1, n_CIA

        ! Check for 3 species special
        if (CIA(s)%i3 .eqv. .True.) then
          ni = 3
        else
          ni = 2
        end if

        do i = 1, ni
          exists = .False.
          do j = 1, n_sp
            if (ni == 2) then
              if (CIA(s)%sp_con(i) == sp_list(j)) then
                CIA(s)%iVMR(i) = j
                exists = .True.
                exit
              end if
            else if (ni == 3) then
              if (CIA(s)%sp_con_3(i) == sp_list(j)) then
                CIA(s)%iVMR_3(i) = j
                exists = .True.
                exit
              end if
            end if
          end do

          if (exists .eqv. .False.) then
            print*, 'ERROR - Specified CIA species component not found in prf VMR list - STOPPING'
            if (ni == 2) then
              print*, 'Species 2 part: ', CIA(s)%sp, CIA(s)%sp_con(i)
            else if (ni == 3) then
              print*, 'Species 3 part: ', CIA(s)%sp, CIA(s)%sp_con_3(i)
            end if
            stop
          end if

        end do
      end do
    end if

    ! Read in prepared CIA tables
    if (n_CIA > 0) then

      do s = 1, n_CIA

        print*, 'Reading CIA species: ', trim(CIA(s)%sp), ' @ ', trim(CIA(s)%path)

        ! Skip H-
        if ((trim(CIA(s)%sp) == 'H-')) then
          cycle
        end if

        if (CIA(s)%form == 1) then

          open(newunit=u,file=trim(CIA(s)%path),action='read',status='old')

          ! Allocate CIA table temperature arrays
          allocate(CIA(s)%T(CIA(s)%nT))

          ! Read and allocate data until error (end of file)
          do n = 1, CIA(s)%nT
            read(u,*,iostat=stat) name, wn_s, wn_e, nrec, temp_r, kmax, dum
            !print*, n, name, wn_s, wn_e, nrec, temp_r, kmax, dum
            if (n == 1) then
              ! Allocate CIA table wn and table value array
              allocate(CIA(s)%wn(nrec))
              allocate(CIA(s)%kap(nrec,CIA(s)%nT))
              CIA(s)%nwl = nrec
            end if
            ! Check if end of file reached
            if (is_iostat_end(stat)) then
              print*,'Reached end of HITRAN CIA file: ', CIA(s)%sp, CIA(s)%path
              exit
            else
              ! Temperature point of table
              CIA(s)%T(n) = temp_r
              ! Read the record data
              do i = 1, nrec
                read(u,*) CIA(s)%wn(i), CIA(s)%kap(i,n)
                !print*, i, CIA(s)%wn(i), CIA(s)%kap(i,n)
              end do
            end if

          end do
          close(u)

        else if (CIA(s)%form == 2) then

          CIA(s)%nT = 8

          if (CIA(s)%sp == 'He-') then
            CIA(s)%nwl = 16
          else if (CIA(s)%sp == 'H2-') then
            CIA(s)%nwl = 18
          end if

          open(newunit=u,file=trim(CIA(s)%path),action='read',status='old')
          read(u,*)

          allocate(CIA(s)%wl(CIA(s)%nwl))
          allocate(CIA(s)%wn(CIA(s)%nwl))
          allocate(CIA(s)%T(CIA(s)%nT))
          allocate(CIA(s)%kap(CIA(s)%nwl,CIA(s)%nT))

          read(u,*) CIA(s)%T(:)

          do n = 1, CIA(s)%nwl
            read(u,*) CIA(s)%wl(n), CIA(s)%kap(n,:)
            !print*, s, n, CIA(s)%wl(n), CIA(s)%kap(n,:)
          end do
          !CIA(s)%wl(:) = CIA(s)%wl(:) * 1.0e-4_dp
          CIA(s)%wn(:) = 1.0_dp/(CIA(s)%wl(:) * 1.0e-8_dp)

          close(u)
        else if (CIA(s)%form == 3) then
           !! Read in reformatted cia table (for memory efficenct models)
            open(newunit=u,file=trim(CIA(s)%path),action='read',status='old')

            read(u,*) CIA(s)%nT, CIA(s)%nwl

            allocate(CIA(s)%wl(CIA(s)%nwl))
            allocate(CIA(s)%wn(CIA(s)%nwl))
            allocate(CIA(s)%T(CIA(s)%nT))
            allocate(CIA(s)%kap(CIA(s)%nwl,CIA(s)%nT))

            read(u,*) CIA(s)%T(:)
            read(u,*) CIA(s)%wn(:)

            read(u,*)

            CIA(s)%wl(:) = 1.0_dp/CIA(s)%wn(:) * 1e-4_dp

            do i = 1, CIA(s)%nT
              read(u,*) CIA(s)%kap(:,i)
            end do

            !! We must reverse the wn array to conform to low to high wavelength
            ! band structure
            allocate(dum_b(CIA(s)%nwl))
            !! Reverse wn and wl arrays
            dum_b(:) = CIA(s)%wn(:)
            do j = 1,  CIA(s)%nwl
              CIA(s)%wn(j) = dum_b(CIA(s)%nwl-j+1)
            end do

            dum_b(:) = CIA(s)%wl(:)
            do j = 1,  CIA(s)%nwl
              CIA(s)%wl(j) = dum_b(CIA(s)%nwl-j+1)
            end do

            do i = 1, CIA(s)%nT
              dum_b(:) = CIA(s)%kap(:,i)
              do j = 1, CIA(s)%nwl
                CIA(s)%kap(j,i) = dum_b(CIA(s)%nwl-j+1)
              end do
            end do

            deallocate(dum_b)
            close(u)           
        end if
      end do
    end if

    !! log10 copies of the CIA tables for the interpolation
    do s = 1, n_CIA
      if (trim(CIA(s)%sp) == 'H-') cycle
      select case (CIA(s)%form)
      case (1)
        CIA(s)%lx = CIA(s)%wn(:)
      case (2)
        CIA(s)%lx = log10(CIA(s)%wl(:))
      end select
      CIA(s)%lT = log10(CIA(s)%T(:))
      CIA(s)%lkap = log10(max(CIA(s)%kap(:,:), 1.0e-99_dp))
      if ((CIA(s)%form /= 3) .and. ((size(CIA(s)%lx) < 3) .or. (size(CIA(s)%lT) < 3))) then
        print*, 'CIA table needs >= 3 points on both axes for the Bezier interpolation: ', trim(CIA(s)%sp)
        stop
      end if
    end do

    if (n_Ray > 0) then
      ! Read Rayleigh namelist
      allocate(Ray(n_Ray))
      allocate(Ray_form(n_Ray),Ray_sp(n_Ray),Ray_paths(n_Ray))
      read(u_nml, nml=Ray_nml)

      Ray(:)%sp = Ray_sp(:)
      Ray(:)%path = Ray_paths(:)

      do s = 1, n_Ray
        exists = .False.
        do j = 1, n_sp
          if (Ray(s)%sp == sp_list(j)) then
            Ray(s)%iVMR = j
            exists = .True.
            exit
          end if
        end do

        if (exists .eqv. .False.) then
          print*, 'ERROR - Specified Rayleigh species component not found in prf VMR list - STOPPING'
          print*, 'Species: ', s, Ray(s)%sp
          stop
        end if
      end do

      ! Read in Rayleigh scattering data
      if (n_Ray > 0) then
        do s = 1, n_Ray
          allocate(Ray(s)%kap(n_b))
          print*, 'Reading Rayleigh data: ', trim(Ray(s)%sp), ' @ ',trim(Ray(s)%path)
          open(newunit=uR,file=trim(Ray(s)%path),action='read',status='old')
          read(uR,*)
          read(uR,*) Ray(s)%kap(:)
          !print*, s, Ray(s)%kap(:)
          close(uR)
        end do
      end if
    end if


    print*, 'Finished reading'

  end subroutine read_gas_data

  ! Adapted CMCRT H- calculation
  subroutine CIA_Hminus(n_sp, s, b, T, Nl, VMR, k_cia)
    implicit none
    ! John (1988) paramaters
    real(dp), dimension(6), parameter :: &
      & An_ff1 = (/518.1021_dp, 472.2636_dp, -482.2089_dp, 115.5291_dp, 0.0_dp, 0.0_dp/), &
      & Bn_ff1 = (/-734.8666_dp, 1443.4137_dp, -737.1616_dp, 169.6374_dp, 0.0_dp, 0.0_dp/), &
      & Cn_ff1 = (/1021.1775_dp, -1977.3395_dp, 1096.8827_dp, -245.6490_dp, 0.0_dp, 0.0_dp/), &
      & Dn_ff1 = (/-479.0721_dp, 922.3575_dp, -521.1341_dp, 114.2430_dp, 0.0_dp, 0.0_dp/), &
      & En_ff1 = (/93.1373_dp, -178.9275_dp, 101.7963_dp, -21.9972_dp, 0.0_dp, 0.0_dp/), &
      & Fn_ff1 = (/-6.4285_dp, 12.3600_dp, -7.0571_dp, 1.5097_dp, 0.0_dp, 0.0_dp/)

    real(dp), dimension(6), parameter :: &
      & An_ff2 = (/0.0_dp, 2483.3460_dp, -3449.8890_dp, 2200.0400_dp, -696.2710_dp, 88.2830_dp/), &
      & Bn_ff2 = (/0.0_dp, 285.8270_dp, -1158.3820_dp, 2427.7190_dp, -1841.4000_dp, 444.5170_dp/), &
      & Cn_ff2 = (/0.0_dp, -2054.2910_dp, 8746.5230_dp, -13651.1050_dp, 8642.9700_dp, -1863.8640_dp/), &
      & Dn_ff2 = (/0.0_dp, 2827.7760_dp, -11485.6320_dp, 16755.5240_dp, -10051.5300_dp, 2095.2880_dp/), &
      & En_ff2 = (/0.0_dp, -1341.5370_dp, 5303.6090_dp, -7510.4940_dp, 4400.0670_dp, -901.7880_dp/), &
      & Fn_ff2 = (/0.0_dp, 208.9520_dp, -812.9390_dp, 1132.7380_dp, -655.0200_dp, 132.9850_dp/)

    real(dp), dimension(6), parameter :: &
      & Cn_bf = (/152.519, 49.534, -118.858, 92.536, -34.194, 4.982 /)

    !real(dp), parameter :: alf = 1.439e8_dp
    real(dp), parameter :: lam_0 = 1.6419_dp

    integer, intent(in) :: s, b, n_sp
    real(dp), dimension(n_sp), intent(in) :: VMR
    real(dp), intent(in) :: T, Nl
    real(dp), intent(out) :: k_cia

    integer :: n
    real(dp) :: kff, kbf, fbf, xbf, sff
    real(dp) :: T5040

    T5040 = 5040.0_dp / T

    ! Do bound free calculation
    if (wl(b) > lam_0) then
      xbf = 0.0_dp
      !kbf = 0.0_dp
    else
      fbf = 0.0_dp
      do n = 1, 6
        fbf = fbf + Cn_bf(n) * (1.0_dp/wl(b) - 1.0_dp/lam_0)**((real(n,dp)-1.0_dp)/2.0_dp)
      end do
      xbf = 1.0e-18_dp * wl(b)**3 * (1.0_dp/wl(b) - 1.0_dp/lam_0)**(3.0_dp/2.0_dp) * fbf
      !kbf = 0.750_dp * T**(-5.0_dp/2.0_dp) * exp(alf/(lam_0*T)) * (1.0_dp - exp(-alf/(wl(b)*T))) * xbf
    end if

    ! Do free free calculation
    sff = 0.0_dp
    if (wl(b) >= 0.3645_dp) then
      do n = 1, 6
        sff = sff + T5040**((real(n,dp)+1.0_dp)/2.0_dp) &
                & * (wl(b)**2*An_ff2(n) + Bn_ff2(n)  + Cn_ff2(n)/wl(b) + Dn_ff2(n)/wl(b)**2 + En_ff2(n)/wl(b)**3 &
                & + Fn_ff2(n)/wl(b)**4)
      end do
      kff = 1.0e-29_dp * sff
    else if ((wl(b) < 0.3645_dp) .and. (wl(b) > 0.1823_dp)) then
      do n = 1, 6
        sff = sff + T5040**((real(n,dp)+1.0_dp)/2.0_dp) &
                & * (wl(b)**2*An_ff1(n) + Bn_ff1(n)  + Cn_ff1(n)/wl(b) + Dn_ff1(n)/wl(b)**2 + En_ff1(n)/wl(b)**3 &
                & + Fn_ff1(n)/wl(b)**4)
      end do
      kff = 1.0e-29_dp * sff
    else
      kff = 0.0_dp
    end if

    ! xbf is in [cm2 molecule-1] and kff is in [cm4 dyne-1] - convert to [cm-1] before CMCRT output
    kbf = xbf * VMR(CIA(s)%iVMR_3(1)) * Nl !! * H- [molecule cm-3]
    !kbf = kbf * (VMR_lay(CIA_tab(s)%iVMR_3(2),z) * N_lay(z) * VMR_lay(CIA_tab(s)%iVMR_3(3),z) * N_lay(z)) * kb * T ! * P(e-) [dyne cm-2] * H [cm-3]
    kff = kff * (VMR(CIA(s)%iVMR_3(2)) * Nl * VMR(CIA(s)%iVMR_3(3)) * Nl) * kb_cgs * T !! * P(e-) [dyne cm-2] * H [cm-3]

    k_cia = kbf + kff

  end subroutine CIA_Hminus

  subroutine CIA_Heminus(n_sp, s, b, T, Nl, VMR, k_Hem)
    implicit none

    integer, intent(in) :: s, b, n_sp
    real(dp), dimension(n_sp), intent(in) :: VMR
    real(dp), intent(in) :: T, Nl
    real(dp), intent(out) :: k_Hem

    real(dp) :: aa, bb, cc, kff

    ! Polynomial coefficents, frequency f in [Hz]
    aa = 3.397e-46_dp + (-5.216e-31_dp + 7.039e-15_dp/freq(b))/freq(b)
    bb = -4.116e-42_dp + (1.067e-26_dp + 8.135e-11_dp/freq(b))/freq(b)
    cc = 5.081e-37_dp + (-8.724e-23_dp - 5.659e-8_dp/freq(b))/freq(b)
    kff = aa*T + bb + cc/T

    kff = kff * VMR(CIA(s)%iVMR(1)) * Nl * VMR(CIA(s)%iVMR(2)) * Nl !  He [cm-3] * e- [cm-3]

    ! End units [cm-1]

    k_Hem = kff

  end subroutine CIA_Heminus

  !! Free-free (form 2) table: Bezier interpolation of log10 kap in log10 wl and log10 theta,
  !! theta = 5040/T, edge values outside the table
  subroutine CIA_ff_Bezier(n_sp, s, b, T, Nl, VMR, kff)
    implicit none

    integer, intent(in) :: s, b, n_sp
    real(dp), dimension(n_sp), intent(in) :: VMR
    real(dp), intent(in) :: T, Nl
    real(dp), intent(out) :: kff

    integer :: iwl, iT
    logical :: ewl, eT
    real(dp) :: lwl, lth, aval

    lwl = log10(wl(b) * 1.0e4_dp)
    lth = log10(5040.0_dp / T)
    call bezier_index(CIA(s)%lx(:), lwl, iwl, ewl)
    call bezier_index(CIA(s)%lT(:), lth, iT, eT)
    aval = 10.0_dp**bezier_2D(CIA(s)%lx(:), CIA(s)%lT(:), CIA(s)%lkap(:,:), lwl, lth, iwl, ewl, iT, eT)

    kff = 1.0e-26_dp * aval * (VMR(CIA(s)%iVMR(1)) * Nl * VMR(CIA(s)%iVMR(2)) * Nl) * kb_cgs * T !! * P(e-) [dyne cm-2] * (H2 or He) [cm-3]

  end subroutine CIA_ff_Bezier

  subroutine CIA_ff(n_sp, s, b, T, Nl, VMR, kff)
    implicit none

    integer, intent(in) :: s, b, n_sp
    real(dp), dimension(n_sp), intent(in) :: VMR
    real(dp), intent(in) :: T, Nl
    real(dp), intent(out) :: kff

    integer :: iwl, iwl1, iT, iT1
    real(dp) :: T5040, wlA
    real(dp) :: xval, x0, x1, yval, y0, y1
    real(dp) :: a00, a10, a01, a11, aval

    !! Outside the table: its edge values (this gave zero below T = 5040/3.6 = 1400 K)
    T5040 = min(max(5040.0_dp / T, CIA(s)%T(1)), CIA(s)%T(CIA(s)%nT))
    wlA = min(max(wl(b) * 1.0e4_dp, CIA(s)%wl(1)), CIA(s)%wl(CIA(s)%nwl))

    call locate(CIA(s)%wl(:),wlA,iwl)
    iwl = min(max(iwl, 1), CIA(s)%nwl - 1)
    iwl1 = iwl + 1

    call locate(CIA(s)%T(:),T5040,iT)
    iT = min(max(iT, 1), CIA(s)%nT - 1)
    iT1 = iT + 1

    !! wn and T are within the table bounds
    xval = wlA ; x0 = CIA(s)%wl(iwl) ; x1 = CIA(s)%wl(iwl1)
    yval = T5040 ; y0 = CIA(s)%T(iT) ; y1 = CIA(s)%T(iT1)
    a00 = CIA(s)%kap(iwl,iT) ; a10 = CIA(s)%kap(iwl1,iT)
    a01 = CIA(s)%kap(iwl,iT1) ; a11 = CIA(s)%kap(iwl1,iT1)

    ! Perform bi-linear interpolation
    call bilinear_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)

   ! print*, s,b,wl(b), T, aval

    kff = 1.0e-26_dp * aval * (VMR(CIA(s)%iVMR(1)) * Nl * VMR(CIA(s)%iVMR(2)) * Nl) * kb_cgs * T !! * P(e-) [dyne cm-2] * (H2 or He) [cm-3]

  end subroutine CIA_ff


  subroutine find_CIA_consituents(n_CIA)
    implicit none

    integer, intent(in) :: n_CIA
    integer :: s

    do s = 1, n_CIA

      select case(CIA(s)%sp)

      case('H2-H2')
        CIA(s)%sp_con(1) = 'H2'
        CIA(s)%sp_con(2) = 'H2'

        CIA(s)%nT = 113

      case('H2-He','He-H2')
        CIA(s)%sp_con(1) = 'H2'
        CIA(s)%sp_con(2) = 'He'

        CIA(s)%nT = 334

      case('H2-H','H-H2')
        CIA(s)%sp_con(1) = 'H2'
        CIA(s)%sp_con(2) = 'H'

        CIA(s)%nT = 4

      case('H-He','He-H')
        CIA(s)%sp_con(1) = 'He'
        CIA(s)%sp_con(2) = 'H'

        CIA(s)%nT = 10

      case('H-')
        CIA(s)%sp_con_3(1) = 'H-'
        CIA(s)%sp_con_3(2) = 'H'
        CIA(s)%sp_con_3(3) = 'e-'

        CIA(s)%i3 = .True.

      case('He-')
        CIA(s)%sp_con(1) = 'He'
        CIA(s)%sp_con(2) = 'e-'

      case('H2-')
        CIA(s)%sp_con(1) = 'H2'
        CIA(s)%sp_con(2) = 'e-'

      case default
        print*, 'ERROR - CIA species constituents could not be found - STOPPING'
        print*, 'Species: ', CIA(s)%sp
        stop
      end select

    end do

  end subroutine find_CIA_consituents

  subroutine locate(arr, var, idx)
    implicit none

    integer, intent(out) :: idx
    real(dp), dimension(:), intent(in) :: arr
    real(dp),intent(in) ::  var
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

  ! Perform linear interpolation  space
  subroutine linear_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    lxval = xval
    lx1 = x1; lx2 = x2
    ly1 = y1; ly2 = y2

    norm = 1.0_dp / (lx2 - lx1)

    yval = ((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)

  end subroutine linear_interp

  subroutine bilinear_interp(xval, yval, x1, x2, y1, y2, a11, a21, a12, a22, aval)
    implicit none

    real(dp), intent(in) :: xval, yval, x1, x2, y1, y2, a11, a21, a12, a22
    real(dp) :: lxval, lyval, lx1, lx2, ly1, ly2, la11, la21, la12, la22
    real(dp), intent(out) :: aval
    real(dp) :: norm

    lxval = xval ; lyval = yval
    lx1 = x1 ; lx2 = x2 ; ly1 = y1 ; ly2 = y2
    la11 = a11 ; la21 = a21 ; la12 = a12 ; la22 = a22

    norm = 1.0_dp / (lx2 - lx1) / (ly2 - ly1)

    aval = la11 * (lx2 - lxval) * (ly2 - lyval) * norm &
      & + la21 * (lxval - lx1) * (ly2 - lyval) * norm &
      & + la12 * (lx2 - lxval) * (lyval - ly1) * norm &
      & + la22 * (lxval - lx1) * (lyval - ly1) * norm

    aval = aval

  end subroutine bilinear_interp

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    lxval = log10(xval)
    lx1 = log10(x1); lx2 = log10(x2)
    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / (lx2 - lx1)

    yval = 10.0_dp**((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)

  end subroutine linear_log_interp

  subroutine bilinear_log_interp(xval, yval, x1, x2, y1, y2, a11, a21, a12, a22, aval)
    implicit none

    real(dp), intent(in) :: xval, yval, x1, x2, y1, y2, a11, a21, a12, a22
    real(dp) :: lxval, lyval, lx1, lx2, ly1, ly2, la11, la21, la12, la22
    real(dp), intent(out) :: aval
    real(dp) :: norm

    lxval = log10(xval) ; lyval = log10(yval)
    lx1 = log10(x1) ; lx2 = log10(x2) ; ly1 = log10(y1) ; ly2 = log10(y2)
    la11 = log10(a11) ; la21 = log10(a21) ; la12 = log10(a12) ; la22 = log10(a22)

    norm = 1.0_dp / (lx2 - lx1) / (ly2 - ly1)

    aval = la11 * (lx2 - lxval) * (ly2 - lyval) * norm &
      & + la21 * (lxval - lx1) * (ly2 - lyval) * norm &
      & + la12 * (lx2 - lxval) * (lyval - ly1) * norm &
      & + la22 * (lxval - lx1) * (lyval - ly1) * norm

    aval = 10.0_dp**(aval)

  end subroutine bilinear_log_interp

  ! Perform Bezier interpolation
  !! ---------------------------------------------------------------------------------------
  !! Shared opacity-table interpolation (identical in Exo-FMS_column_ck and Exo-FMS_2).
  !! Monotonic quadratic Bezier interpolation on 3-point stencils, applied in log10 space
  !! by the callers (log10 k against log10 T, log10 p, ...). Outside a table, or on its end
  !! points, the edge value is used; a table axis with fewer than 3 points is interpolated
  !! linearly (1D only - 2D tables must have >= 3 points on both axes).
  !! ---------------------------------------------------------------------------------------

  !! Stencil start i1 of xg(i1:i1+2) for x with xg(1) < x < xg(n) (edge = .false.), or the
  !! edge index i1 = 1 or n (edge = .true.) when x is on or outside the end points.
  !! Interior: xg(i1+1) < x <= xg(i1+2), except the first interval, xg(1) < x <= xg(2).
  subroutine bezier_index(xg, x, i1, edge)
    implicit none

    real(dp), dimension(:), intent(in) :: xg
    real(dp), intent(in) :: x
    integer, intent(out) :: i1
    logical, intent(out) :: edge

    integer :: n

    n = size(xg)
    edge = .True.
    if (x <= xg(1)) then
      i1 = 1
    else if (x >= xg(n)) then
      i1 = n
    else
      edge = .False.
      call locate(xg, x, i1)
      i1 = min(max(i1 - 1, 1), n - 2)
    end if

  end subroutine bezier_index

  !! Value at x of the 1D table f(xg): Bezier inside (linear if size(xg) < 3), edge values outside
  real(dp) function bezier_1D(xg, f, x) result(val)
    implicit none

    real(dp), dimension(:), intent(in) :: xg, f
    real(dp), intent(in) :: x

    integer :: i1
    logical :: edge

    if (size(xg) == 1) then
      val = f(1)
    else if (size(xg) == 2) then
      if (x <= xg(1)) then
        val = f(1)
      else if (x >= xg(2)) then
        val = f(2)
      else
        val = f(1) + (f(2) - f(1)) * (x - xg(1)) / (xg(2) - xg(1))
      end if
    else
      call bezier_index(xg, x, i1, edge)
      if (edge) then
        val = f(i1)
      else
        call Bezier_interp(xg(i1:i1+2), f(i1:i1+2), 3, x, val)
      end if
    end if

  end function bezier_1D

  !! Value at (x, y) of the 2D table f(xg, yg), with the stencils (ix, ex) and (iy, ey) from
  !! bezier_index: Bezier in y at the three x stencil points, then in x; edge values outside
  real(dp) function bezier_2D(xg, yg, f, x, y, ix, ex, iy, ey) result(val)
    implicit none

    real(dp), dimension(:), intent(in) :: xg, yg
    real(dp), dimension(:,:), intent(in) :: f
    real(dp), intent(in) :: x, y
    integer, intent(in) :: ix, iy
    logical, intent(in) :: ex, ey

    integer :: i
    real(dp), dimension(3) :: fx

    if (ex .and. ey) then
      val = f(ix,iy)
    else if (ex) then
      call Bezier_interp(yg(iy:iy+2), f(ix,iy:iy+2), 3, y, val)
    else if (ey) then
      call Bezier_interp(xg(ix:ix+2), f(ix:ix+2,iy), 3, x, val)
    else
      do i = 1, 3
        call Bezier_interp(yg(iy:iy+2), f(ix+i-1,iy:iy+2), 3, y, fx(i))
      end do
      call Bezier_interp(xg(ix:ix+2), fx(:), 3, x, val)
    end if

  end function bezier_2D

  !! Monotonic quadratic Bezier interpolation through 3 points: the control point is set
  !! from the weighted 3-point derivative at xi(2); if that control point is not strictly
  !! between the interval's end values (overshoot, a local extremum or a flat interval),
  !! the one-sided slope is used instead, which gives linear interpolation. x outside
  !! [xi(1), xi(3)] gives the end value (no extrapolation).
  subroutine Bezier_interp(xi, yi, ni, x, y)
    implicit none

    integer, intent(in) :: ni
    real(dp), dimension(ni), intent(in) :: xi, yi
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y

    real(dp) :: dx, dx1, dy, dy1, ww, yc, t

    dx = xi(2) - xi(1)
    dx1 = xi(3) - xi(2)
    dy = yi(2) - yi(1)
    dy1 = yi(3) - yi(2)

    if (x <= xi(2)) then
      ! left hand side interpolation
      ww = dx1/(dx + dx1)
      yc = yi(2) - dx/2.0_dp * (ww*dy/dx + (1.0_dp - ww)*dy1/dx1)
      if ((yc <= min(yi(1),yi(2))) .or. (yc >= max(yi(1),yi(2)))) then
        yc = yi(2) - dx/2.0_dp * (dy/dx)
      end if
      t = min(max((x - xi(1))/dx, 0.0_dp), 1.0_dp)
      y = (1.0_dp - t)**2 * yi(1) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(2)
    else
      ! right hand side interpolation
      ww = dx/(dx + dx1)
      yc = yi(2) + dx1/2.0_dp * (ww*dy1/dx1 + (1.0_dp - ww)*dy/dx)
      if ((yc <= min(yi(2),yi(3))) .or. (yc >= max(yi(2),yi(3)))) then
        yc = yi(2) + dx1/2.0_dp * (dy1/dx1)
      end if
      t = min(max((x - xi(2))/dx1, 0.0_dp), 1.0_dp)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine Bezier_interp

  subroutine sort2(N,RA,RB)
    integer, intent(in) :: N
    integer :: L, IR, I, J
    real(dp), dimension(N), intent(inout) :: RA, RB
    real(dp) :: RRA, RRB
    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
      L=L-1
      RRA=RA(L)
      RRB=RB(L)
    ELSE
      RRA=RA(IR)
      RRB=RB(IR)
      RA(IR)=RA(1)
      RB(IR)=RB(1)
      IR=IR-1
      IF(IR.EQ.1)THEN
        RA(1)=RRA
        RB(1)=RRB
        RETURN
      ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
      IF(J.LT.IR)THEN
        IF(RA(J).LT.RA(J+1)) THEN
          J=J+1
        ENDIF
      ENDIF
      IF(RRA.LT.RA(J))THEN
        RA(I)=RA(J)
        RB(I)=RB(J)
        I=J
        J=J+J
      ELSE
        J=IR+1
      ENDIF
      GO TO 20
    ENDIF
    RA(I)=RRA
    RB(I)=RRB
    GO TO 10
  end subroutine sort2

  real(dp) function BB(wl_in, T_in)
    implicit none

    ! Planck function in wavelength units

    real(dp), intent(in) :: wl_in, T_in
    real(dp) :: left, right
    real(dp) :: wl_cm

    wl_cm = wl_in * 1.0e-4_dp

    left = (2.0_dp * hp * c_s**2)/(wl_cm)**5
    right = 1.0_dp / (exp((hp * c_s) / (wl_cm * kb * T_in)) - 1.0_dp)
    BB = left * right

  end function BB

end module ck_opacity_mod
