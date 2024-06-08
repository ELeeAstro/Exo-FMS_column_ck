import numpy as np
import matplotlib.pylab as plt

data = np.loadtxt('FMS_RC_flx.out')
Pe = data[:,1]
sw_up = data[:,2]
sw_down = data[:,3]
sw_net = data[:,4]
lw_up = data[:,5]
lw_down = data[:,6]
lw_net = data[:,7]

fig = plt.figure()

plt.plot(sw_up,Pe/1e5,ls='dashed',lw=3,label='sw_up')
plt.plot(sw_down,Pe/1e5,lw=3,ls='dashed',label='sw_down')
plt.plot(sw_net,Pe/1e5,ls='dashed',lw=3,label='sw_net')

plt.plot(lw_up,Pe/1e5,lw=3,label='lw_up')
plt.plot(lw_down,Pe/1e5,lw=3,label='lw_down')
plt.plot(lw_net,Pe/1e5,lw=3,label='lw_net')

plt.ylabel(r'Pressure [bar]')
plt.xlabel(r'Flux [W m$^{-2}$]')
plt.legend()

plt.yscale('log')
plt.gca().invert_yaxis()

yticks = [1000,100,10,1,0.1,0.01,1e-3,1e-4,1e-5,1e-6]
yticks_lab = ['1000','100','10','1','0.1','0.01','10$^{-3}$','10$^{-4}$','10$^{-5}$','10$^{-6}$']

plt.legend()

plt.ylim(1000,1e-6)
plt.yticks(yticks,yticks_lab)


plt.show()

