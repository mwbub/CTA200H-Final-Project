import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia

# Note: convert this program into a function.

# Stars to query
stars = ["GJ 440", "Lacaille 8760", "Vega"]

# Get data on each star from SIMBAD, for use in epoch adjustment
customSimbad = Simbad()
customSimbad.add_votable_fields('pmra', 'pmdec', 'plx', 'rv_value')
simbad_table = customSimbad.query_objects(stars)

# Convert ra and dec from SIMBAD into SkyCoord objects for each star
simbad_coords = [SkyCoord(simbad_table['RA'][i], simbad_table['DEC'][i], 
                   unit = (u.hourangle, u.deg)) for i in range(len(stars))]

# Hold the data on each star in a tuple of the form 
# (ra, dec, parallax, pmra, pmdec, radial_velocity)
simbad_vals = [(simbad_coords[i].ra.value, simbad_coords[i].dec.value, 
                simbad_table['PLX_VALUE'][i], simbad_table['PMRA'][i], 
                simbad_table['PMDEC'][i], simbad_table['RV_VALUE'][i])
                for i in range(len(stars))]

# Choose a circle radius to query in Gaia
radius = (1.0*u.arcsec).to(u.deg)

# Perform a cone search for each star with epoch adjustment
queries = ["""
SELECT *
FROM gaiadr2.gaia_source
WHERE 1=CONTAINS(
        POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),
        CIRCLE('ICRS', 
               COORD1(EPOCH_PROP_POS({0},{1},{2},{3},{4},{5},2000,2015.5)),
               COORD2(EPOCH_PROP_POS({0},{1},{2},{3},{4},{5},2000,2015.5)),
               {6}))
""".format(*simbad_vals[i], radius.value) for i in range(len(stars))]

# Perform each query and get the result in a table
jobs = [Gaia.launch_job_async(query) for query in queries]
results = [job.get_results() for job in jobs]
