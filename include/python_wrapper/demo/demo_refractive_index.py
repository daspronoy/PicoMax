# %%
from picomax import *
import matplotlib.pyplot as plt


# %%
c = picomax()
c.crystal = 'zincblende'
c.a = 4.35
c.vg = np.array([3,4,8,11,-0.225,0.005,0.118,0.127,-0.545,-0.005,0.118,0.121])
c.nvalence = 4
c.refpoint = 'G'
c.qpt = np.array([[0,0,0]])
c.nband = 20
c.nomega = 1001
c.domega = 0.02
c.omega2 = 2
c.epsilon = 0.3
c.kk_transform = 1
c.kptopt = 1
c.kptorder = 15
c.nchi = 1
c.encut = 200
c.wdir = '/home/mjh92'

c.run(switch='Xij')

# %%

plt.figure(figsize=(4,3.5))
ax = plt.axes()
plt.xlim([0, 40])
line = plt.plot(c.omega,1+np.real(c.Xij[:,0,0,0,0]),'b')
plt.setp(line,linewidth=0.75)
line = plt.plot(c.omega,np.imag(c.Xij[:,0,0,0,0]),'r')
plt.setp(line,linewidth=0.75)
plt.xlabel('Frequency (eV)')
plt.ylabel('Permittivity (a.u.)')

# %%
