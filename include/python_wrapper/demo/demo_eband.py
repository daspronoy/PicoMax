
# %%
from picomax import *
import matplotlib.pyplot as plt


# %%
c = picomax()
c.set(crystal='zincblende',
      a=4.35,
      vg=np.array([3,4,8,11,-0.225,0.005,0.118,0.127,-0.545,-0.005,0.118,0.121]),
      nvalence=4,
      refpoint='G',
      wdir='/home/mjh92/',
      qpath='L,G,X,K,G',
      qnum=[20,20,10,20],
      nband=30)

# run EPM calculation
c.run(switch='eband')





# %%

plt.figure(figsize=(4,3.5))
ax = plt.axes()
for n in range(c.nband):
    line = plt.plot(c.iqpt,c.eband[:,n],'b')
    plt.setp(line,linewidth=0.75)

plt.ylabel('Energy (eV)')
plt.xlim([0, c.nqpt-1])
plt.ylim([-20, 20])
plt.grid(True, axis='x', color='k', alpha=0.5, linestyle='--')
ax.set_xticks(np.append([0],np.cumsum(c.qnum)),minor=False)
ax.set_xticklabels(c.qpath.split(','))
ax.xaxis.grid(True,which='major')
# %%


# plt.savefig('/home/mjh92/test.svg')