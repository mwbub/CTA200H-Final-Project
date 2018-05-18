import os
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import gaia_tools.load as load
import gaia_tools.xmatch as xmatch
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014

# Load the APOGEE red clump catalogue
apogeerc = load.apogeerc()

print("Loading data...")

# Cross match APOGEE red clump with Gaia DR2
if not os.path.exists("apogeerc_xmatched.csv"):
    xmatched = xmatch.cds(apogeerc, xcat = "vizier:I/345/gaia2", 
                          savefilename = "apogeerc_xmatched.csv")[0]
else:
    xmatched = xmatch.cds_load("apogeerc_xmatched.csv")

# Filter out stars with no radial velocity
rv_mask = [xmatched['radial_velocity'] != -9999.99]
xmatched_rv = xmatched[rv_mask]

ra = xmatched_rv['ra']
dec = xmatched_rv['dec']
d = 1/xmatched_rv['parallax']
pmra = xmatched_rv['pmra']
pmdec = xmatched_rv['pmdec']
rv = xmatched_rv['radial_velocity']

nsteps = 100
physical = False
t = np.linspace(0, 20, nsteps)
coord_path = "integrated_{}_steps".format(nsteps)

if physical:
    t *= u.Gyr
    coord_path += "_physical"
    
coord_path += ".npy"

if not os.path.exists(coord_path):
    print("Generating orbits...")
    orbits = [Orbit(vxvv = [ra[i], dec[i], d[i], pmra[i], pmdec[i], rv[i]], 
                radec = True) for i in range(len(xmatched_rv))]
    
    print("Integrating orbits...")
    for orbit in orbits:
        orbit.integrate(t, MWPotential2014)
        
    print("Getting coordinates...")
    x = [[orbit.x(time) for orbit in orbits] for time in t]
    y = [[orbit.y(time) for orbit in orbits] for time in t]
    np.save(coord_path, (x,y))
else:
    x, y = np.load(coord_path)

print("Plotting...")
plt.scatter(x[-1], y[-1], s = 0.1)
plt.xlim(-30,30)
plt.ylim(-30,30)
plt.show()