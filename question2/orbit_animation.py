import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

nsteps = 1000
physical = True
filepath = "data/integrated_{}_steps".format(nsteps)
if physical:
    filepath += "_physical"
filepath += ".npy"

if not os.path.exists('animations'):
    os.mkdir('animations')

if os.path.exists(filepath):
    x, y = np.load(filepath)
else:
    raise OSError
    
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set(xlim = (-30, 30), ylim = (-30, 30), xlabel = "x (kpc)", 
       ylabel = "y (kpc)", title = "APOGEE Red Clump Integrated Over 20 Gyr")
scatter = ax.scatter(x[0], y[0], s=0.1)   

def animate(i):
    scatter.set_offsets(np.c_[x[i], y[i]])

anim = animation.FuncAnimation(fig, animate, interval=100, frames=len(x)-1)
plt.draw()
plt.show()

anim.save('animations/' + os.path.splitext(os.path.basename(filepath))[0] + 
          "_animated.gif", writer='imagemagick')
