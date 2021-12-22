# Exo-FMS_column_nongrey

Elspeth KH Lee - May 2021

This is one part of a series of codes that build upon different two-stream approaches and schemes, primarily useful for the GCM modeling community.
This is the non-grey, picket fence version, with three bands in the visible, representing incident radiation from a host star, and two bands in the IR, representing the internal radiation propagating inside the planetary atmosphere.

Some useful references useful for non-grey, picket fence modeling are: \
Chandrasekhar (1935) \
Parmentier & Guillot (2014) \
Parmentier et al.  2015) \
Lee et al. (2021)

The non-grey, picket fence scheme builds upon the semi-grey scheme, by having 3 bands in the visible, and 2 bands in the IR.
The 3 V bands attempt to better represent the deposition of stellar flux into the atmosphere, and therefore produce a more realistic stellar heating profile.
The 2 IR bands (picket fence scheme), are used to both represent the line opacity in one band, and the continuum opacity in the other.
The fraction of the total flux in each band is controlled by a parameter (Beta_IR).
This allows the atmosphere to cool from two photospheric regions, rather than one in the semi-grey scheme, leading to more realistic cooling in the upper atmosphere.
However, at very low pressures, isothermal T-p profiles are still produced.

Using the fitting functions and relations from Parmentier et al. (2014, 2015), a quite realistic T-p profile can be produced, since the fitting functions were designed to reproduce the results of a correlated-k scheme as best as possible.
Unfortunately the numerical results in 1D can in certain circumstances produce small spiky profiles, primarily due to the sensitivity of the opacity fitting function to the temperature, or, when the fitting function is outside it's valid range of temperatures (~100-4000 K) or pressures (~1e-4-1e3 pa).
Inside a GCM these small spikes are typically smoothed out by the dynamical processes.

This code performs various two-stream approaches from the literature in a non-grey, picket fence context:
1. Isothermal layer approximation
2. Toon et al. method
3. Short Characteristics method
4. Heng et al. method
5. Neil Lewis' scattering code, following Pierrehumbert (2010)
6. Mendonca et al. method (IN DEVELOPMENT)
7. Two-stream DISORT version (w. modifications by Xianyu Tan)

This emulates a single column inside the Exo-FMS GCM and is useful for testing and developing new techniques
as they would perform inside a GCM setting. This is also useful to see differences in each method and their various approximations.

We also include a dry convective adjustment schemes, currently only 'Ray_adj', based on Raymond Pierrehumbert's python code.

# Namelist options

In the file 'FMS_RC.nml' you can select different options that control the simulation

ts_scheme: \
'Isothermal' - Isothermal ts method \
'Isothermal_2' - Isothermal ts method - high optical depth version \
'Toon' - Toon et al. ts method \
'Toon_scatter' - Toon et al. ts method with scattering (IN DEVELOPMENT) \
'Shortchar' -  Short characteristics method \
'Heng' - Heng et al. method \
'Lewis_scatter' - Neil Lewis's scattering code, following Pierrehumbert (2010) \
'Lewis_scatter_sw' - 'Lewis' but only shortwave scattering ('Shortchar' for IR) \
'Mendonca' - Mendonca et al. method (IN DEVELOPMENT) \
'Disort_scatter' - two-stream DISORT version with scattering

opac_scheme: \
'Freedman' - Uses the Rosseland mean fitting function from Freedman et al. (2014) \
'Valencia' - Uses the Rosseland mean fitting function from Valencia et al. (2013) \

adj_scheme: \
'Ray_dry' - Ray Pierrehumbert's dry convective adjustment scheme

The option 'None' for each of these scheme will it off (e.g. To run without conv adjustment set adj_scheme = 'None')

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
AB - Bond albedo \
fl - The Heng et al. (2011) parameter used for pressure dependent IR optical depths

sw_ac(3) - shortwave single scattering albedo (constant at all layers) \
sw_gc(3) - shortwave asymmetry factor (constant at all layers) \
lw_ac(2) - longwave single scattering albedo (constant at all layers) \
lw_gc(2) - longwave asymmetry factor (constant at all layers)

iIC - Initial condition selection integer (4 = Parmentier et al. (2015) profile using Tirr and Tint) \
corr - Flag to perform the adiabatic gradient correction in the initial conditions \
met - metallicty in dex solar (M/H) \
table_num - 1 = Parmentier et al. (2015) table with TiO & VO, 2 = without TiO & VO table

# Plotting results

A python plotting routine, 'plot_TP.py' is provided to plot the results of the simulation.

# Gaussian Ordinance values

In methods that use Gaussian quadrature to perform the mu integration in intensity (to calculate the flux), various values and weights can be changed for testing at the top of the two-stream modules.
You will need to clean and recompile the code if these are changed.

# Personal recommendations

For non-scattering problems, we generally recommend that the short characteristics method be used, as it is fast, efficient, very stable and also very accurate. This is currently what is used inside Exo-FMS for the Hot Jupiter simulations, and is even fast enough for high-resolution cases.
For scattering problems we recommend the two stream DISORT version, it is very reliable but generally slower compared to other scattering methods.

# Future developments

We will include versions that include multiple-scattering in the future. \
Additional opacity schemes from the literature can be added. \
Ability to include solid surface temperatures and temperature evolution, this involves some extra switches and boundary conditions. \
Window fractions and functions.
Interpolating directly from the Freedman et al. (2014) tables rather than using the fitting function.
