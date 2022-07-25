import numpy as np

#Common physical constants - CODATA 2020
kb = 1.380649e-16 # erg K-1 - Boltzmann's constant
hp = 6.62607015e-27 # erg s - Planck's constant
c_s = 2.99792458e10 # cm s-1 - Vacuum speed of light
m_el = 9.1093837015e-28 # g - electron mass
a_fine = 7.2973525693e-3 # fine structure constant
Comp_e = hp/(m_el*c_s) # Compton wavelength for electrons
sigT = ((8.0*np.pi)/3.0) * ((a_fine * Comp_e)/(2.0*np.pi))**2 # Thomson cross section [cm2]

#Parameters for H
wl_ly = 121.567 * 1.0e-7 # Lyman alpha wavelength [cm]
f_ly = c_s/wl_ly
w_l = (2.0 * np.pi * f_ly) / 0.75
cp = [1.26537,3.73766,8.8127,19.1515, \
    39.919,81.1018,161.896,319.001,622.229,1203.82]


# Code that generates Rayleigh xsections for various species
wl = np.loadtxt('../wavelengths_GCM_32.txt',skiprows=1)
nwl = len(wl) - 1
wl = (wl[:1] + wl[1:])/2.0
wl_A = wl * 1e4
wl_cm = wl / 1e4
wn = 1.0/wl_cm
freq = c_s*wn

H_xsec = np.zeros(nwl)
# Calculate for H
for l in range(nwl):
  w = 2.0 * np.pi * freq[l]
  wwl = w/w_l

  #Lee and Kim (2004)
  if (wwl <= 0.6):
    #Low energy limit
    xsec = 0.0
    for p in range(10):
      xsec = xsec + (cp[p] * wwl**(2 * p))
      xsec = xsec * wwl**4
  else:
    #High energy limit (approaching Lyman alpha wavelengths)
     wb = (w - 0.75*w_l)/(0.75*w_l)
     xsec = (0.0433056/wb**2)*(1.0 - 1.792*wb - 23.637*wb**2 - 83.1393*wb**3 \
       - 244.1453*wb**4 - 699.473*wb**5)
  #Multiply by Thomson x-section
  H_xsec[l] = xsec * sigT
  #print(l,wl[l],H_xsec[l])

# Calculate for H2
# Use Irwin (2009)+ paramaters (same as Cox 2000)
H2_xsec = np.zeros(nwl)
for l in range(nwl):
  A = 13.58e-5
  B = 7.52e-3
  n_ref = A * (1.0 + B/wl[l]**2) + 1.0
  King = 1.0 
  nd_stp = 2.65163e19
  H2_xsec[l] = ((24.0 * np.pi**3 * wn[l]**4)/(nd_stp**2)) \
             * ((n_ref**2 - 1.0)/(n_ref**2 + 2.0))**2  * King
  #print(l,wl[l],H2_xsec[l])

# Calculate for He
# Use Thalman et al. (2014) parameters
He_xsec = np.zeros(nwl)
for l in range(nwl):
  A = 2283.0
  B = 1.8102e13
  C = 1.5342e10
  n_ref = (A + (B / (C - wn[l]**2)))/1.0e8 + 1.0
  King = 1.0
  nd_stp = 2.546899e19
  He_xsec[l] = ((24.0 * np.pi**3 * wn[l]**4)/(nd_stp**2)) \
             * ((n_ref**2 - 1.0)/(n_ref**2 + 2.0))**2  * King
  #print(l,wl[l],He_xsec[l])

# Calculate for e-
e_xsec = np.zeros(nwl)
e_xsec[:] = sigT

# Output H
output_name = 'Ray_H_'+str(nwl)+'.txt'
output = open(output_name,'w')
output.write('H' + '\n')
output.write(' '.join(str(g) for g in H_xsec[:]) + '\n')
output.close()

# Output H2
output_name = 'Ray_H2_'+str(nwl)+'.txt'
output = open(output_name,'w')
output.write('H2' + '\n')
output.write(' '.join(str(g) for g in H2_xsec[:]) + '\n')
output.close()

# Output He
output_name = 'Ray_He_'+str(nwl)+'.txt'
output = open(output_name,'w')
output.write('He' + '\n')
output.write(' '.join(str(g) for g in He_xsec[:]) + '\n')
output.close()

# Output e-
output_name = 'Ray_e-_'+str(nwl)+'.txt'
output = open(output_name,'w')
output.write('e-' + '\n')
output.write(' '.join(str(g) for g in e_xsec[:]) + '\n')
output.close()

