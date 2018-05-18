import os
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import gaia_tools.load as load
import gaia_tools.xmatch as xmatch
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014
from galpy.util.bovy_conversion import time_in_Gyr

# Load the APOGEE red clump catalogue
apogeerc = load.apogeerc()

print("Loading data...")

if not os.path.exists("data"):
    os.mkdir("data")

# Cross match APOGEE red clump with Gaia DR2
if not os.path.exists("data/apogeerc_xmatched.csv"):
    xmatched = xmatch.cds(apogeerc, xcat = "vizier:I/345/gaia2", 
                          savefilename = "data/apogeerc_xmatched.csv")[0]
else:
    xmatched = xmatch.cds_load("data/apogeerc_xmatched.csv")

# Filter out stars with no radial velocity
rv_mask = (xmatched['radial_velocity'] != -9999.99)
xmatched_rv = xmatched[rv_mask]

nstars = len(xmatched_rv)
ra = xmatched_rv['ra']
dec = xmatched_rv['dec']
d = 1/xmatched_rv['parallax']
pmra = xmatched_rv['pmra']
pmdec = xmatched_rv['pmdec']
rv = xmatched_rv['radial_velocity']

nsteps = 1000
t = np.linspace(0, 20, nsteps)
coord_path = "data/integrated_{}_steps".format(nsteps)

physical = False
if physical:
    t *= u.Gyr
    coord_path += "_physical"
else:
    t /= time_in_Gyr(220., 8.)
    
coord_path += ".npy"

if not os.path.exists(coord_path):
    print("Generating orbits...")
    orbits = [Orbit(vxvv = [ra[i], dec[i], d[i], pmra[i], pmdec[i], rv[i]], 
                    ro = 8., vo = 220., radec = True) for i in range(nstars)]
    
    print("Integrating orbits...")
    for orbit in orbits:
        orbit.integrate(t, MWPotential2014)
        
    print("Getting coordinates...")
    x = np.stack([orbit.x(t) for orbit in orbits], axis = 1)
    y = np.stack([orbit.y(t) for orbit in orbits], axis = 1)
    np.save(coord_path, (x,y))
else:
    x, y = np.load(coord_path)

print("Plotting...")
x_coords = x.flatten()
y_coords = y.flatten()
H, xedges, yedges = np.histogram2d(x_coords, y_coords, bins = 500, 
                                   range = ((-15, 15), (-15, 15)))
X, Y = np.meshgrid(xedges, yedges)
plt.pcolormesh(X, Y, H, cmap = 'plasma')
plt.show()
