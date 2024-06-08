# Exo-FMS_column_ck

[HEAVY RECONSTRUCTION IN PROGRESS - TODO add references and descriptions and add comments]

Elspeth KH Lee - Feb 2022

Major Update History:
 - Jun 2022 - overhaul of opacities + Bezier methods
 - Jun 2024 - major overhaul of code and new methods added

This is the third part of a series of codes that builds upon different two-stream approaches and schemes, primarily useful for the GCM modelling community.
This is the correlated-k (corr-k) version, currently set up for 11, 30 or 32 bands as discussed in Kataria et al. (2013), Showman et al. (2009) and Amundsen et al. (2014) respectively for modelling hot Jupiter atmospheres.

Some useful references useful for corr-k modelling are: \
Goody et al. (1989) \
Lacis & Oinas (1991) \
Showman et al. (2009) \
Kataria et al. (2013) \
Grimm & Heng (2015) \
Amundsen et al. (2014, 2017) \
Lee et al. (2021) \
and plenty of others...

Unfortunately the numerical results in 1D can in certain circumstances produce small spiky profiles, primarily due to issues either relating to the opacity table fidelity, important optical absorbing species condensing (at low pressure e.g. TiO) leading to large opacity gradients, interpolation error or vertical grid resolution.
However, inside a GCM these spikes are typically smoothed out by the dynamical processes.

To compile enter 'make' in the main directory. To remove compiled code enter 'make clean'.
Some compiler options for gfortran, nvfortran and ifort are provided in the makefile.

This code performs various atmospheric radiative-transfer approaches from the literature using corr-k, split into short and long wave methods.

For the shortwave method we have:
1. Direct beam only method
2. Adding method
3. Spherical harmonic adding method (SDA)
4. Toon89 multiple scattering method
5. Spherical harmonic two-stream multiple scattering method (SH2)
6. Spherical harmonic four-stream multiple scattering method (SH4)
7. Two-stream DISORT version

For the longwave method we have:
1. Absorption approximation with exponential in tau function (AA_E)
2. Absorption approximation with linear in tau function (AA_L)
3. Short characteristics (sc) with linear interpolants 
4. Toon89 multiple scattering method
5. Variational iteration method (VIM)

This emulates a single column inside the Exo-FMS GCM and is useful for testing and developing new techniques
as they would perform inside a GCM setting. This is also useful to see differences in each method and their various approximations.

For longwave radiation integrate the Planck function we utilize the method described here (Widger, W. K. and Woodall, M. P. 1976): \
https://spectralcalc.com/blackbody/inband_radiance.html \
We found this to be very quick and stable method, producing similar results to the Disort version but with quicker convergence in most cases, especially with the way it is coded here by using the results of the previous band to get the next bin value. \
An alternative would be to produce a table of integrated fluxes for a temperature range and for each wavelength band and interpolate to that, as is done in other GCM models e.g. SPARC/MITgcm.

For the longwave radiation, we use the weighted essentially non oscillatory order 4 interpolation method (WENO4) from Janett et al. (2019). This allows for highly smooth interpolation of the temperature from the layers to the levels (typically this is required in GCMs) and strong gradients in temperature are accurately interpolated through. Linear extrapolation is used for the end points to avoid numerical issues and overshoot.

We also include dry convective adjustment schemes, currently 'Ray_adj', based on Raymond Pierrehumbert's python code and mixing length theory (MLT).

# Requirements & Opacity data:

You need to download CIA tables from the HITRAN database, mainly H2-H2, H2-He, H2-H and He-H for gas giants. \
Typically you add them to the 'cia' data folder.

Our pre-mixed k-tables contained in this repository have over 30 species of interest (including UV-OPT species, Fe, Fe+, SiO, TiO and VO) to exoplanet atmospheres and cool dwarf stars. \
The pre-mixed tables are valid between a temperature of 100-6100 K and pressure 1e-8 - 1000 bar, making them suitable for UHJ modelling too. \
Importantly we include the important species responsible for the upper atmosphere temperature inversion as found in Lothringer et al. (2020). \
Alternative pre-mixed tables without UV-OPT absorbers labeled 'nUVOPT' in the ck data directory. \
The pre-mixed tables were calculated assuming local condensation only (i.e. species that have condensed are at their saturation point (S=1) VMRs)

More individual k-tables can be produced upon request and will be made available publicly at some point, right now we include some important species for 11 bins (in the directory "11_ck_data_g8" since it is much faster) for the random overlap resort rebin (RORR) and adaptive equivalent extinction (AEE) methods.

# Namelist options

### See FMS_RC.nml_RORR for a random overlap resort rebin example using Burrows analytic CE abundances.
### See FMS_RC.nml_AEE for a adaptive equivalent extinction example using Burrows analytic CE abundances.
### See FMS_RC.nml_PM for a pre-mixed example using interpolation from the CE grid.

In the file 'FMS_RC.nml' you can select different options that control the simulation, a default set-up for HD 209458b is provided.

### &FMS_RC_nml

sw_scheme: \
'sw_direct' - direct beam only \ 
'sw_adding' - adding method \
'sw_SDA' - SDA method \
'sw_Toon' - Toon89 method \
'sw_SH2' - SH2 method \
'sw_SH4' - SH4 method \ 
'sw_disort_ts' - two stream disort method

lw_scheme: \
'lw_AA_E' - AA_E method \
'lw_AA_L' - AA_L method \
'lw_sc_linear' - sc linear method \
'lw_VIM' - VIM method \
'lw_Toon' - Toon89 method \ 
'lw_disort_ts' - two steam disort method

opac_scheme: \
'ck' - Use the corr-k scheme here (only option)

adj_scheme: \
'Ray_dry' - Ray Pierrehumbert's dry convective adjustment scheme

CE_scheme: \
'interp' - interpolation from CE table \

The option 'None' for each of these scheme will it off (e.g. To run without conv adjustment set adj_scheme = 'None')

nb - number of wavelength bands \
wl_sh - file with wavelength edges \
stellarf_sh - file with stellar surface fluxes

Rs - Stellar radius (Solar units) \
sm_ax - semi-major axis (AU)

n_ck - number of corr-k species (Pre-mixed = 1)
ng  - number of g-ordinance values in k-tables \
nsp - number of species that require VMRs \
n_cia - number of CIA species \
n_ray - number of Rayligh scattering species

nlay, a_sh , b_sh - the number of layers, and filenames that contain the a and b constants for the hybrid sigma grid

pref - reference surface pressure (pa)

t_step - time step in seconds \
nstep - number of integer timesteps \
Rd_air - specific gas constant of the air (J kg-1 K-1)\
cp_air - specific heat capacity (constant pressure) of the air (J kg-1 K-1) \
grav - gravitational acceleration constant (m s-2) \
mu_z - cosine angle of the solar zenith angle \
Tirr - Irradiation temperature \
Tint - Internal temperature

k_V - visible band opacity (m2 kg-1) \
k_IR - IR band opacity (m2 kg-1) \
fl - The Heng et al. (2011) parameter used for pressure dependent IR optical depths \
met - metallicty in dex solar (M/H)

zcorr - include zenith angle correction (.False.) \
zcorr_meth - zenith angle correction method (1,2)  \
radius - radius of the planet at surface (m)

iIC - Initial condition selection integer (4 = Parmentier et al. (2015) profile using Tirr and Tint) \
corr - Flag to perform the adiabatic gradient correction in the initial conditions \
table_num - 1 = Parmentier et al. (2015) table with TiO & VO, 2 = without TiO & VO table

### &sp_nml

sp_list - list of species that need a VMR \
VMR_tab_sh - path to VMR tables

NOTE: Specifc ordering of species must be used for Burrows and the interpolation tables (See examples)

### &gw_nml

gw - the Gaussian weights for the k-tables

### &ck_nml

PM - .True. (for pre-mixed tables), .False. \
RORR - .True. (for random overlap resort rebin), .False. \
AEE - .True. (for adaptive equivalent extinction), .False. \
ck_interp -  1 = bilinear interpolate ck tables, 2 = Bezier interpolate ck tables \
ck_form  - format of k-tables (2 = Helios-k format used here) \
ck_sp - names of species for corr-k routines ('PM' for pre-mixed) \
ck_paths - path to the corr-k tables in order of species in ck_sp

### &cia_nml

cia_form - format for CIA table (1 = HITRAN, 2 = He-, H2- tables, 0 = special e.g. H- opacity) \
cia_sp - names of CIA species \
cia_paths - paths to the CIA tables in order of cia_sp, use '.' for special species

### &ray_nml

ray_form - format of Rayleigh scattering table (only = 1 currently) \
ray_sp - names of Rayleigh scattering species \
ray_paths - paths to Rayleigh scattering tables in order of ray_sp

# Plotting results

A python plotting routine, 'plot_TP.py' is provided to plot the results of the simulation.

# Gaussian Ordinance values

In methods that use Gaussian quadrature to perform the mu integration in intensity (to calculate the flux), various values and weights can be changed for testing at the top of the two-stream modules.
You will need to clean and recompile the code if these are changed. Beware, some codes use the classical Gauss-Legendre method (flux = 2pi * intensity), while others include the weighting in the intensity (flux = pi * intensity).

# Personal recommendations

For shortwave scattering problems we recommend the adding method as included, or if more accuracy is needed the SDA or Toon89 methods.
For longwave scattering problems we recommend the linear absorption approximation method or if more accuracy is required the VIM or Toon89 methods.


# Future developments

Ability to include solid surface temperatures and temperature evolution, this involves some extra switches and boundary conditions. \
Improve interpolation of opacity tables to reduce spikiness. \
Add general cloud opacity features. \ 
