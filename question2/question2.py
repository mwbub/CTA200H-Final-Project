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

if not os.path.exists("data/apogeerc_xmatched.csv"):
    # Cross match APOGEE red clump with Gaia DR2
    xmatched = xmatch.cds(apogeerc, xcat = "vizier:I/345/gaia2", 
                          savefilename = "data/apogeerc_xmatched.csv")[0]
else:
    # Load the cross matched table if already computed
    xmatched = xmatch.cds_load("data/apogeerc_xmatched.csv")

# Filter out stars with no radial velocity
rv_mask = (xmatched['radial_velocity'] != -9999.99)
xmatched_rv = xmatched[rv_mask]

# Separate the columns of xmatched_rv
ra = xmatched_rv['ra']
dec = xmatched_rv['dec']
d = 1/xmatched_rv['parallax']
pmra = xmatched_rv['pmra']
pmdec = xmatched_rv['pmdec']
rv = xmatched_rv['radial_velocity']

# Choose the number of time steps for the orbit integration
nsteps = 1000
t = np.linspace(0, 20, nsteps)
coord_path = "data/integrated_{}_steps".format(nsteps)

# Choose whether to use astropy physical units, or to convert time
# to Gyr using galpy.util.bovy_conversion
physical = False
if physical:
    t *= u.Gyr
    coord_path += "_physical"
else:
    t /= time_in_Gyr(220., 8.)
    
coord_path += ".npy"

if not os.path.exists(coord_path):
    # Generate and integrate the orbits for each star
    print("Generating orbits...")
    orbits = [Orbit(vxvv = [ra[i], dec[i], d[i], pmra[i], pmdec[i], 
                            rv[i]], ro = 8., vo = 220.,
                            radec = True) for i in range(len(xmatched_rv))]
    
    print("Integrating orbits...")
    for orbit in orbits:
        orbit.integrate(t, MWPotential2014)
    
    # Collect and save the x and y coordinates for each star at each time step
    print("Getting coordinates...")
    x = np.stack([orbit.x(t) for orbit in orbits], axis = 1)
    y = np.stack([orbit.y(t) for orbit in orbits], axis = 1)
    np.save(coord_path, (x,y))
else:
    # Load the orbit coordinates if already computed
    x, y = np.load(coord_path)

# Make a density plot of the orbit coordinates
print("Plotting...")
x_coords = x.flatten()
y_coords = y.flatten()
H, xedges, yedges = np.histogram2d(x_coords, y_coords, bins = 250,
                                   range = ((-15, 15), (-15, 15)))
X, Y = np.meshgrid(xedges, yedges)
plt.pcolormesh(X, Y, H, cmap = 'plasma')
plt.show()
