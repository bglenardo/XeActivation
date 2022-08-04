import numpy as np
from matplotlib import pyplot as plt


flux = 1e4 # n/cm^2/s, value from MNRC
xs = 3e-24 # cm^2, estimated from Xe126 capture xs at ~a few x 10^-8 eV.
Xe127_tau = 52.4 * 24*60*60 # time constant in seconds

irrad_time = 200 # days

natXe_mass = 1 # kg

natXe_126frac = 0.00089

natXe_amu = 131.2

natXe_molperkg = 1/(natXe_amu * 0.001)

Xe126_atmperkg = 6.02e23 * natXe_126frac * natXe_molperkg
print('Xe126 atoms per kg of Xe: {}',Xe126_atmperkg)

# Array containing the absolute number of Xe127 atoms in a
# kg of natXe in 1 second intervals
secs = irrad_time * 24 * 60 * 60
N127 = np.zeros(500000)
time = np.linspace(0.,secs,500000)
dt = time[1] - time[0]

for i in range(1,len(N127)):
  N127[i] = N127[i-1] +\
              (-N127[i-1]/Xe127_tau +\
              (Xe126_atmperkg - N127[i-1]) * xs * flux -\
              N127[i-1] * xs * flux) * dt
  #print(' {} '.format(i))
  #print('Term 1: {}'.format(-N127[i-1]/Xe127_tau))
  #print('Term 2: {}'.format((Xe126_atmperkg - N127[i-1]) * xs * flux))
  #print('Term 3: {}'.format(N127[i-1] * xs * flux))
  #print('N127: {}'.format(N127[i]))

plt.plot(time/60/60/24,N127,'-r')
plt.xscale('linear')
plt.yscale('log')
plt.xlabel('Time (days)')
plt.ylabel('Number of Xe127 atoms in sample')
plt.show() 
