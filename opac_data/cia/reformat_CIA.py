import numpy as np

sp = ['H2-H2','H2-He','H2-H','He-H']
nsp = len(sp)

# Load in H2-H2 CIA
fname = 'H2-H2_2011.cia'
f = open(fname,'r')
nT_H2 = 113
T_H2 = np.zeros(nT_H2)
kap_H2 = np.zeros((nT_H2,9981))
wn_H2 = np.arange(20,10001,1)
for i in range(nT_H2):
  h = f.readline().split()
  wn_s = float(h[1])
  wn_e = float(h[2])
  nrec = int(h[3])
  T_H2[i] = float(h[4])
  #print(i,T_H2[i])
  for j in range(nrec):
      h = f.readline().split()
      kap_H2[i,j] = np.log10(float(h[1]))

# Load in H2-He CIA
fname = 'H2-He_2011.cia'
f = open(fname,'r')
nT_He = 334
T_He = np.zeros(nT_He)
kap_He = np.zeros((nT_He,19981))
wn_He = np.arange(20,20001,1)
for i in range(nT_He):
  h = f.readline().split()
  wn_s = float(h[1])
  wn_e = float(h[2])
  nrec = int(h[3])
  T_He[i] = float(h[4])
  #print(i,T_He[i])
  for j in range(nrec):
      h = f.readline().split()
      kap_He[i,j] = np.log10(float(h[1]))

# Load in H2-H CIA
fname = 'H2-H_2011.cia'
f = open(fname,'r')
nT_H = 4 
T_H = np.zeros(nT_H)
kap_H = np.zeros((nT_H,9901))
wn_H = np.arange(100,10001,1)
for i in range(nT_H):
  h = f.readline().split()
  wn_s = float(h[1])
  wn_e = float(h[2])
  nrec = int(h[3])
  T_H[i] = float(h[4])
  #print(i,T_He[i])
  for j in range(nrec):
      h = f.readline().split()
      kap_H[i,j] = np.log10(float(h[1]))

# Load in He-H CIA
fname = 'He-H_2011.cia'
f = open(fname,'r')
nT_HeH = 10
T_HeH = np.zeros(nT_HeH)
kap_HeH = np.zeros((nT_HeH,10951))
wn_HeH = np.arange(50,11001,1)
for i in range(nT_HeH):
  h = f.readline().split()
  wn_s = float(h[1])
  wn_e = float(h[2])
  nrec = int(h[3])
  T_HeH[i] = float(h[4])
  #print(i,T_He[i])
  for j in range(nrec):
      h = f.readline().split()
      kap_HeH[i,j] = np.log10(float(h[1]))

# Read in wavenumber grid edges
wl = np.loadtxt('../wavelengths_GCM_32.txt',skiprows=1)
nwl = len(wl)
nb = len(wl) - 1


wnb = 1.0/(wl / 1e4)
wnb = wnb[::-1]
wnb_c = np.zeros(nb)
wnb_c = (wnb[0:-1] + wnb[1:])/2.0


print(wnb)
print(wnb_c)


k_wn_H2H2 = np.zeros((nT_H2,nb))
k_wn_H2He = np.zeros((nT_He,nb))
k_wn_H2H = np.zeros((nT_H,nb))
k_wn_HeH = np.zeros((nT_HeH,nb))

for s in range(nsp):
  print(s,sp[s])
  if (sp[s] == 'H2-H2'):
    for t in range(nT_H2):
      k_wn_H2H2[t,:] = 10.0**np.interp(wnb_c[:],wn_H2[:],kap_H2[t,:],left=-99,right=-99)
  elif (sp[s] == 'H2-He'):
    for t in range(nT_He):
      k_wn_H2He[t,:] = 10.0**np.interp(wnb_c[:],wn_He[:],kap_He[t,:],left=-99,right=-99)
  elif (sp[s] == 'H2-H'):
    for t in range(nT_H):
      k_wn_H2H[t,:] = 10.0**np.interp(wnb_c[:],wn_H[:],kap_H[t,:],left=-99,right=-99)
  elif (sp[s] == 'He-H'):
    for t in range(nT_HeH):
      k_wn_HeH[t,:] = 10.0**np.interp(wnb_c[:],wn_HeH[:],kap_HeH[t,:],left=-99,right=-99)

print('Outputting')

fname_H2 = 'H2-H2_reform_'+str(nb)+'.txt'
f = open(fname_H2,'w')
f.write(str(nT_H2) + ' ' +  str(nb) + ' ' + '\n')
f.write(" ".join(str(g) for g in T_H2[:]) + '\n')
f.write(" ".join(str(g) for g in wnb_c[:]) + '\n')
f.write('\n')
for t in range(nT_H2):
    f.write(" ".join(str(g) for g in k_wn_H2H2[t,:]) + '\n')
f.close()

fname_He = 'H2-He_reform_'+str(nb)+'.txt'
f = open(fname_He,'w')
f.write(str(nT_He) + ' ' +  str(nb) + ' ' + '\n')
f.write(" ".join(str(g) for g in T_He[:]) + '\n')
f.write(" ".join(str(g) for g in wnb_c[:]) + '\n')
f.write('\n')
for t in range(nT_He):
    f.write(" ".join(str(g) for g in k_wn_H2He[t,:]) + '\n')
f.close()

fname_H = 'H2-H_reform_'+str(nb)+'.txt'
f = open(fname_H,'w')
f.write(str(nT_H) + ' ' +  str(nb) + ' ' + '\n')
f.write(" ".join(str(g) for g in T_H[:]) + '\n')
f.write(" ".join(str(g) for g in wnb_c[:]) + '\n')
f.write('\n')
for t in range(nT_H):
    f.write(" ".join(str(g) for g in k_wn_H2H[t,:]) + '\n')
f.close()

fname_HeH = 'He-H_reform_'+str(nb)+'.txt'
f = open(fname_HeH,'w')
f.write(str(nT_HeH) + ' ' +  str(nb) + ' ' + '\n')
f.write(" ".join(str(g) for g in T_HeH[:]) + '\n')
f.write(" ".join(str(g) for g in wnb_c[:]) + '\n')
f.write('\n')
for t in range(nT_HeH):
    f.write(" ".join(str(g) for g in k_wn_HeH[t,:]) + '\n')
f.close()
