import numpy as np
import matplotlib.pylab as plt

# E. K.H. Lee - Nov 2020

# A simple code to generate a hybrid sigma grid for FMS - add the coefficents to fv_eta.F90 and modify the arrays in the fortan code to use 
# Used to generate the 53 220 - 1e-6 bar HJ grid and 1000 - 1e-6 bar deep HJ grid
# p0 = 2.2e7, pu = 5e-2, ps = 1e1, nlev = 53

# Generally, the idea is that the linear spacing at lower pressure is more mutable to radiative transfer calculations.
# The optical depth spacing will be larger, so you do not have many cells requiring small tau approximations.
# In addition, exending the grid to lower pressure (e.g. 1e-6) is useful for RT post-processing (e.g. Transmission spectra)
# as well as post-processing for photo-chemical modelling (e.g. VULCAN).

# The linear part will adjust the hybrid grid to higher pressure than a pure sigma. 
# So pu generally has to be at lower pressure than the target lowest pressure
# For example 5e-2 rather than 1e-1 for a top level near 1e-1 pa
# Zoom in on the scatter plot at lowest pressure to check if the target pressure is reached and adjust as appropriate

p0 = 1.0e7 # 0 altitude level pressure [pa]
pu = 1e1  # Highest altitude level pressure [pa]
ps = 9e0   # Approx pressure of linear spacing start [pa] 

# Number of levels
nlev = 41
nlay = nlev - 1

# Log-space pressure from p0 to pu
p = np.logspace(np.log10(pu),np.log10(p0),nlev,endpoint=True)

# Linear-space pressure from ps to pu
a = np.linspace(pu,ps,nlev)

# Last a has to be exactly 0
a[-1] = 0.0

# Last sigma point has to be exactly 1
p[-1] = p0

# Calculate coefficents and also find pressure levels for plotting
sigma = np.zeros(nlev) 
pure_sigma = np.zeros(nlev)
hyb_sigma = np.zeros(nlev)
idx = np.zeros(nlev)
for i in range(nlev):
    idx[i] = i + 1
    sigma[i] = p[i]/p0
    pure_sigma[i] = sigma[i]*p0
    hyb_sigma[i] = a[i] + sigma[i]*p0


# Write a and b coefficents to a file
fname = 'Hybrid_sigma_grid_a_'+str(nlay)+'.txt'
f = open(fname,'w')

for i in range(nlev):
  #f.write(str(a[i]) +  ',' + '\n')
  f.write(str(a[i]) + '\n')

fname = 'Hybrid_sigma_grid_b_'+str(nlay)+'.txt'
f = open(fname,'w')

for i in range(nlev):
  #f.write(str(sigma[i]) +  ',' + '\n')
  f.write(str(sigma[i])  + '\n')



## Figure plotting - plot both pure sigma and hybrid sigma

fig = plt.figure()

plt.scatter(idx,pure_sigma/1e5,c='black',label='Pure')
plt.scatter(idx,hyb_sigma/1e5,c='red',label='Hybrid')

plt.legend()

plt.yscale('log')

plt.gca().invert_yaxis()

plt.show()
