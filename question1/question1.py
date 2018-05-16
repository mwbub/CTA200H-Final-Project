import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia

def star_to_orbit(star_name, radius = 1.0*u.arcsec):
    # Get data on the star from SIMBAD, for use in epoch adjustment
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields('pmra', 'pmdec', 'plx', 'rv_value')
    simbad_table = custom_simbad.query_object(star_name)
    
    # Convert ra and dec from SIMBAD into a SkyCoord object
    simbad_coord = SkyCoord(simbad_table['RA'][0], 
                            simbad_table['DEC'][0], 
                            unit = (u.hourangle, u.deg))
    
    # Store the SIMBAD values in a tuple of the form 
    # (ra, dec, plx, pmra, pmdec, rv)
    simbad_vals = (simbad_coord.ra.value,
                   simbad_coord.dec.value,
                   simbad_table['PLX_VALUE'][0],
                   simbad_table['PMRA'][0],
                   simbad_table['PMDEC'][0],
                   simbad_table['RV_VALUE'][0])
    
    # Circle radius to query Gaia, in degrees
    radius_deg = radius.to(u.deg)
    
    # Perform a cone search with epoch adjustment
    query =  """
             SELECT source_id, ra, dec, pmra, pmdec, parallax, radial_velocity
             FROM gaiadr2.gaia_source
             WHERE 1=CONTAINS(
             POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),
             CIRCLE('ICRS', 
             COORD1(EPOCH_PROP_POS({0},{1},{2},{3},{4},{5},2000,2015.5)),
             COORD2(EPOCH_PROP_POS({0},{1},{2},{3},{4},{5},2000,2015.5)), {6}))
             """.format(*simbad_vals, radius_deg.value)
             
    # Perfrom the cone search and get the result in a table
    job = Gaia.launch_job_async(query)
    return job.get_results()

