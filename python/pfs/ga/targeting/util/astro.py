import numpy as np

def ABmag_to_flux(mag, zero=0):
    """Convert AB mags to flux"""
    return 10**(-0.4 * (mag - zero + 48.6))

def flux_to_ABmag(flux, zero=0):
    """Convert flux to AB mags"""
    return -2.5 * np.log10(flux) - 48.6 + zero

def ABmag_to_sigma(mag, conversion, sky, sky_sigma=None, zero=0, softening=0):
    if sky_sigma is None:
        sky_sigma = sky
    
    flux = ABmag_to_flux(mag, zero=zero)
    counts = flux * conversion
    rel_error = np.sqrt(counts + sky + sky_sigma ** 2) / counts + softening
    mag_error = 2.5 * rel_error / np.log(10)
    return mag_error

def ABmag_sigma_to_conversion(mag, mag_error, sky, sky_sigma=None, zero=0, softening=0):
    """
    Estimate the conversion function from a measures magnitude and quoted error.
    """

    if sky_sigma is None:
        sky_sigma = sky
    
    flux = ABmag_to_flux(mag, zero=zero)
    rel_error = mag_error * np.log(10) / 2.5 - softening
    counts = 1 / rel_error ** 2 + np.sqrt(1 / rel_error ** 4 + 4 * (sky - sky_sigma ** 2) / rel_error **2)
    conversion = counts / flux
    return conversion

def sdss_to_hsc(sdss_g, sdss_r, sdss_i, sdss_z):
    hsc_g = sdss_g - 0.00816446 - 0.08366937 * (sdss_g - sdss_r) - 0.00726883 * (sdss_g - sdss_r)**2
    hsc_i = sdss_i + 0.00130204 - 0.16922042 * (sdss_i - sdss_z) - 0.01374245 * (sdss_i - sdss_z)**2
    return hsc_g, hsc_i

def dm_to_kpc(dm):
    return 1e-3 * 10**((dm + 5.0) / 5.0)

def kpc_to_dm(kpc):
    return 5 * np.log(kpc * 1e3) / np.log(10) - 5

def arcmin_to_kpc(s, r):
    return 