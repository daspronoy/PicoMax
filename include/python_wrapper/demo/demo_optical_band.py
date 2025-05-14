# %%
from picomax import *
import matplotlib.pyplot as plt


# %%
c = picomax()
c.crystal = 'zincblende'
c.a = 4.35
c.vg = np.array([3,4,8,11,-0.225,0.005,0.118,0.127,-0.545,-0.005,0.118,0.121])
c.nvalence = 4
c.refpoint = np.array([0,0,0])
c.qpt = np.array([[0,0,0]])
c.nband = 20
c.nomega = 2
c.domega = 0.02
c.omega2 = 2
c.epsilon = 0.2
c.kk_transform = 0
c.kptopt = 1
c.kptorder = 11
c.nchi = 9
c.encut = 200
c.wdir = '/home/mjh92'

c.run(switch='Xij')

# %%
