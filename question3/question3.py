import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from galpy.potential import MWPotential2014
from galpy.actionAngle import actionAngleStaeckel

def icrs_to_EccZmaxRperiRap(ra: tuple, dec: tuple, d: tuple, pmra: tuple, 
                            pmdec: tuple, rv: tuple) -> tuple:
    """ Converts ICRS coordinates and their Gaussian uncertainties for a star
    in the Milky War into eccentricity, maximum height above the plane, 
    pericentre, and apocentre with uncertainties.
    
    Parameters:
        
        ra - RA with uncertainty, given as a tuple of the form (RA, RA_ERR).
        Can be Quantity, otherwise given in degrees.
         
        dec - Dec with uncertainty, given as a tuple of the form 
        (DEC, DEC_ERR). Can be Quantity, otherwise given in degrees.
        
        d - Distance with uncertainty, given as a tuple of the form (D, D_ERR).
        Can be Quantity, otherwise given in kpc.
        
        pmra - Proper motion in RA with uncertainty, given as a tuple of the 
        form (PMRA, PMRA_ERR). Can be Quantity, otherwise given in mas/yr.
    
        pmdec - Proper motion in DEC with uncertainty, given as a tuple of the 
        form (PMDEC, PMDEC_ERR). Can be Quantity, otherwise given in mas/yr.
    
        rv - Radial velocity with uncertainty, given as a tuple of the form
        (RV, RV_ERR). Can be Quantity, otherwise given in km/s.
    
    Returns:
        
        (e, zmax, rperi, rap), each a tuple of the form (VALUE, VALUE_ERR).
    """
    # Convert each value and its error into astropy Quantities
    ra_val, ra_err = u.Quantity(ra, u.deg)
    dec_val, dec_err = u.Quantity(dec, u.deg)
    d_val, d_err = u.Quantity(d, u.kpc)
    pmra_val, pmra_err = u.Quantity(pmra, u.mas/u.yr)
    pmdec_val, pmdec_err = u.Quantity(pmdec, u.mas/u.yr)
    rv_val, rv_err = u.Quantity(rv, u.km/u.s)
    
    # Get parameters for a 6D normal distribution to sample
    mean = np.array([ra_val.value, dec_val.value, d_val.value, pmra_val.value, 
                     pmdec_val.value, rv_val.value])
    cov = np.diag([ra_err.value, dec_err.value, d_err.value, pmra_err.value, 
                    pmdec_err.value, rv_err.value])
    
    # Sample the normal distribution
    samples = np.random.multivariate_normal(mean, cov, size=100000)
    ra_samples = samples[:, 0] * u.deg
    dec_samples = samples[:, 1] * u.deg
    d_samples = samples[:, 2] * u.kpc
    pmra_samples = samples[:, 3] * (u.mas/u.yr)
    pmdec_samples = samples[:, 4] * (u.mas/u.yr)
    rv_samples = samples[:, 5] * (u.km/u.s)
    
    coords = SkyCoord(frame="icrs", ra=ra_samples, dec=dec_samples,
                      distance=d_samples, pm_ra_cosdec=pmra_samples, 
                      pm_dec=pmdec_samples, radial_velocity = rv_samples)
    
    # Convert to galactocentric frame
    gal_coords = coords.transform_to("galactocentric")
    gal_coords.representation_type = "cylindrical"
    
    # Convert each value into consistent units
    R = gal_coords.rho.to(u.kpc)
    phi = gal_coords.phi.to(u.deg)
    z = gal_coords.z.to(u.kpc)
    vR = gal_coords.d_rho.to(u.km/u.s)
    vT = (gal_coords.d_phi * R).to(u.km/u.s, 
         equivalencies = u.dimensionless_angles())
    vz = gal_coords.d_z.to(u.km/u.s)
    
    # Calculate EccZmaxRperiRap
    aAS = actionAngleStaeckel(pot=MWPotential2014, delta=0.4)
    orbit_vals = aAS.EccZmaxRperiRap(R, vR, vT, z, vz, phi)
    
    # Get the mean and standard deviation of the samples
    ecc_val, zmax_val, rperi_val, rap_val = np.mean(orbit_vals, axis=1)
    ecc_err, zmax_err, rperi_err, rap_err = np.std(orbit_vals, axis=1)
    
    return ((ecc_val, ecc_err),
            (zmax_val, zmax_err),
            (rperi_val, rperi_err),
            (rap_val, rap_err))
    
# Apply the function to an example star from Gaia
if __name__ == '__main__':
    import warnings
    warnings.filterwarnings("ignore")
    
    job = Gaia.launch_job_async("""
            SELECT TOP 1
            ra, ra_error, dec, dec_error, parallax, 
            parallax_error, pmra, pmra_error, pmdec,
            pmdec_error, radial_velocity, radial_velocity_error
            FROM gaiadr2.gaia_source 
            WHERE radial_velocity IS NOT NULL
            AND parallax_error < parallax
            AND ABS(radial_velocity) < 100.0
            """)
    stars = job.get_results()
    
    ra = (stars[0]["ra"], stars[0]["ra_error"])
    dec = (stars[0]["dec"], stars[0]["dec_error"])
    pmra = (stars[0]["pmra"], stars[0]["pmra_error"])
    pmdec = (stars[0]["pmdec"], stars[0]["pmdec_error"])
    rv = (stars[0]["radial_velocity"], stars[0]["radial_velocity_error"])
    
    d_val = 1/stars[0]["parallax"]
    d_err = abs((stars[0]["parallax_error"]/stars[0]["parallax"]) * d_val)
    d = (d_val, d_err)
    
    coords=icrs_to_EccZmaxRperiRap(ra, dec, d, pmra, pmdec, rv)
    for line in coords:
        print(line)

    
    