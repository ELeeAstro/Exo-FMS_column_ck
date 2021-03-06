# Exo-FMS_column_ck

Elspeth KH Lee - Feb 2022

Major Update History:
 - Jun 2022 - overhaul of opacities + Bezier methods

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

This code performs various two-stream approaches from the literature in a non-grey, picket fence context:
1. Toon et al. method (w. scattering version)
2. Short Characteristics method (w. linear or Bezier interpolants)
3. Two-stream DISORT version (w. modifications by Xianyu Tan)

This emulates a single column inside the Exo-FMS GCM and is useful for testing and developing new techniques
as they would perform inside a GCM setting. This is also useful to see differences in each method and their various approximations.

For the shortwave fluxes, for methods that do not contain a shortwave scattering mode we include the 'adding method' (Mendonca et al. 2015 + references). We detect if any albedo is present in the column, and perform the adding method to calculate the scattered flux, otherwise if there is no albedo only the direct beam is used.

To integrate the Planck function we utilise the method described here (Widger, W. K. and Woodall, M. P. 1976): \
https://spectralcalc.com/blackbody/inband_radiance.html \
We found this to be very quick and stable method, producing similar results to the Disort version but with quicker convergence in most cases, especially with the way it is coded here by using the results of the previous band to get the next bin value. \
An alternative would be to produce a table of integrated fluxes for a temperature range and for each wavelength band and interpolate to that, as is done in other GCM models e.g. SPARC/MITgcm.

We also include a dry convective adjustment schemes, currently only 'Ray_adj', based on Raymond Pierrehumbert's python code.

# Requirements & Opacity data:

You need to download CIA tables from the HITRAN database, mainly H2-H2, H2-He, H2-H and He-H for gas giants. \
Typically you add them to the 'cia' data folder.

Our pre-mixed k-tables contained in this repository have over 30 species of interest (including UV-OPT species, Fe, Fe+, SiO, TiO and VO) to exoplanet atmospheres and cool dwarf stars. \
The pre-mixed tables are valid between a temperature of 200-6100 K and pressure 1e-8 - 1000 bar, making them suitable for UHJ modelling too. \
Importantly we include the important species responsible for the upper atmosphere temperature inversion as found in Lothringer et al. (2020). \
Alternative pre-mixed tables without UV-OPT absorbers labeled 'nUVOPT' in the ck data directory. \
The pre-mixed tables were calculated assuming local rainout (i.e. species that have condensed are at their saturation point (S=1) VMRs)

More individual k-tables can be produced upon request and will be made available publicly at some point, right now we include some important species for 11 bins (in the directory "11_ck_data_g8" since it is much faster) for the random overlap resort rebin (RORR) and adaptive equivalent extinction (AEE) methods.

# Namelist options

### See FMS_RC.nml_RORR for a random overlap resort rebin example using Burrows analytic CE abundances.
### See FMS_RC.nml_AEE for a adaptive equivalent extinction example using Burrows analytic CE abundances.
### See FMS_RC.nml_PM for a pre-mixed example using interpolation from the CE grid.

In the file 'FMS_RC.nml' you can select different options that control the simulation, a default set-up for HD 209458b is provided.

### &FMS_RC_nml

ts_scheme: \
'Toon_scatter' - Toon et al. ts method with scattering \
'Shortchar_linear' -  Short characteristics method with linear interpolants \
'Disort_scatter' - two-stream DISORT version with scattering

opac_scheme: \
'ck' - Use the corr-k scheme here (only option)

adj_scheme: \
'Ray_dry' - Ray Pierrehumbert's dry convective adjustment scheme

CE_scheme: \
'Burrows' - use Burrows analytic CE abundances \
'Interp_bilinear' - use bilinear interpolation from CE table \
'Interp_Bezier' - use Bezier interpolation from CE table \
'Min' - use CE minimisation scheme (in development)

The option 'None' for each of these scheme will it off (e.g. To run without conv adjustment set adj_scheme = 'None')

data_dir - path to data directory \
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

Bezier - use Bezier interpolation for Temperature layer to level interpolation (.True.)

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
You will need to clean and recompile the code if these are changed.

# Personal recommendations

For non-scattering problems, we generally recommend that the short characteristics method be used (either linear or Bezier interpolants), as it is fast, efficient, very stable and also very accurate. This is currently what is used inside Exo-FMS for the Hot Jupiter simulations.
For shortwave scattering problems we recommend the adding method as included (or using the two-stream Toon or DISORT methods), the adding method is generally fast and accurate (enough), especially for grey opacity problems.
For longwave scattering problems we recommend the two stream Toon scattering version.
If that fails try the DISORT code, it is very reliable but generally slower compared to other scattering methods.

# Future developments

We will include versions that include multiple-scattering in the longwave in the future. \
Ability to include solid surface temperatures and temperature evolution, this involves some extra switches and boundary conditions.
