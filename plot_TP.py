import numpy as np
import matplotlib.pylab as plt

data1 = np.loadtxt('FMS_RC_ic.out')
Pic = data1[:,1]
Tic = data1[:,2]

data2 = np.loadtxt('FMS_RC_pp.out')
Ppp = data2[:,1]
Tpp = data2[:,2]
dT_rad = data2[:,3]
dT_conv = data2[:,4]
Kzz = data2[:,5]

fig = plt.figure()

plt.plot(Tic,Pic/1e5,ls='dashed',lw=3,label='Initial Conditions',c='orange')
plt.plot(Tpp,Ppp/1e5,lw=3,label='Numerical Result',c='blue')

plt.ylabel('Pressure [bar]')
plt.xlabel('Temperature [K]')
plt.legend()

plt.yscale('log')
plt.gca().invert_yaxis()

yticks = [1000,100,10,1,0.1,0.01,1e-3,1e-4,1e-5,1e-6]
yticks_lab = ['1000','100','10','1','0.1','0.01','10$^{-3}$','10$^{-4}$','10$^{-5}$','10$^{-6}$']

plt.legend()

plt.ylim(1000,1e-6)
plt.yticks(yticks,yticks_lab)


fig = plt.figure()

plt.plot(dT_rad,Ppp/1e5,lw=3,label='dT_rad',c='red')
plt.plot(dT_conv,Ppp/1e5,lw=3,label='dT_conv',c='blue')
plt.ylabel('Pressure [bar]')
plt.xlabel('dT [K s$^{-1}$]')

plt.yscale('log')
plt.gca().invert_yaxis()

yticks = [1000,100,10,1,0.1,0.01,1e-3,1e-4,1e-5,1e-6]
yticks_lab = ['1000','100','10','1','0.1','0.01','10$^{-3}$','10$^{-4}$','10$^{-5}$','10$^{-6}$']

plt.legend()

plt.ylim(1000,1e-6)
plt.yticks(yticks,yticks_lab)


fig = plt.figure()

plt.plot(Kzz*1e4,Ppp/1e5,lw=3,label='Kzz',c='red')
plt.ylabel(r'Pressure [bar]')
plt.xlabel(r'K$_{\rm zz}$ [cm$^{2}$ s$^{-1}$]')

plt.yscale('log')
plt.xscale('log')
plt.gca().invert_yaxis()

yticks = [1000,100,10,1,0.1,0.01,1e-3,1e-4,1e-5,1e-6]
yticks_lab = ['1000','100','10','1','0.1','0.01','10$^{-3}$','10$^{-4}$','10$^{-5}$','10$^{-6}$']

plt.legend()

plt.ylim(1000,1e-6)
plt.yticks(yticks,yticks_lab)


plt.show()
