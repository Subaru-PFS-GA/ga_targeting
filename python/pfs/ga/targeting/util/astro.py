import numpy as np

EMPTY = object()

def ABmag_to_flux(mag, mag_err=EMPTY, zero=0):
    """Convert AB mags to flux in cgs units"""
    flux = 10**(-0.4 * (mag - zero + 48.6))
    if mag_err is not EMPTY:
        if mag_err is not None:
            flux_err = 0.4 * np.log(10) * flux * mag_err
            return flux, flux_err
        else:
            return flux, None
    else:
        return flux

def ABmag_to_Jy(mag, mag_err=EMPTY, zero=0):
    """Convert AB mags to flux in nJy"""
    flux = 10**(-0.4 * (mag - zero - 8.90))
    if mag_err is not EMPTY:
        if mag_err is not None:
            flux_err = 0.4 * np.log(10) * flux * mag_err
            return flux, flux_err
        else:
            return flux, None
    else:
        return flux

def ABmag_to_nJy(mag, mag_err=EMPTY, zero=0):
    """Convert AB mags to flux in nJy"""
    flux = 10**(-0.4 * (mag - zero - 8.90) + 9)
    if mag_err is not EMPTY:
        if mag_err is not None:
            flux_err = 0.4 * np.log(10) * flux * mag_err
            return flux, flux_err
        else:
            return flux, None
    else:
        return flux

def flux_to_ABmag(flux, flux_err=EMPTY, zero=0):
    """Convert flux in cgs units to AB mags"""
    mag = -2.5 * np.log10(flux) - 48.6 + zero
    if flux_err is not EMPTY:
        if flux_err is not None:
            mag_err = 2.5 / np.log(10) * flux_err / flux
            return mag, mag_err
        else:
            return mag, None
    else:
        return mag

def Jy_to_ABmag(flux, flux_err=EMPTY, zero=0):
    """Convert flux in nJy to AB mags"""
    mag = -2.5 * np.log10(flux) + 8.90 + zero
    if flux_err is not EMPTY:
        if flux_err is not None:
            mag_err = 2.5 / np.log(10) * flux_err / flux
            return mag, mag_err
        else:
            return mag, None
    else:
        return mag

def nJy_to_ABmag(flux, flux_err=EMPTY, zero=0):
    """Convert flux in nJy to AB mags"""
    mag = -2.5 * np.log10(flux * 1e-9) + 8.90 + zero
    if flux_err is not EMPTY:
        if flux_err is not None:
            mag_err = 2.5 / np.log(10) * flux_err / flux
            return mag, mag_err
        else:
            return mag, None
    else:
        return mag

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
    """
    Convert SDSS photometry to HSC using formulae from Komiyama et al. ApJ 853 29K (2018)
    """

    # https://hsc.mtk.nao.ac.jp/pipedoc/pipedoc_7/colorterms.html
    # https://hsc.mtk.nao.ac.jp/pipedoc/pipedoc_8/colorterms.html
    # https://python.hotexamples.com/site/file?hash=0x84cbf64609a6d47c06c3209b08b3186662df8606a83b9cab81665be727fabc57

    # hscPipe 4 colorterms
    hsc_g = sdss_g - 0.00816446 - 0.08366937 * (sdss_g - sdss_r) - 0.00726883 * (sdss_g - sdss_r)**2
    hsc_r = sdss_r + 0.00231810 + 0.01284177 * (sdss_r - sdss_i) - 0.03068248 * (sdss_r - sdss_i)**2
    hsc_i = sdss_i + 0.00130204 - 0.16922042 * (sdss_i - sdss_z) - 0.01374245 * (sdss_i - sdss_z)**2

    # hscPipe 8 colorterms
    # hsc_g = sdss_g - 0.009777 - 0.077235 * (sdss_g - sdss_r) - 0.013121 * (sdss_g - sdss_r)**2
    # hsc_r = sdss_r - 0.000711 - 0.006847 * (sdss_r - sdss_i) - 0.035110 * (sdss_r - sdss_i)**2
    # hsc_i = sdss_i + 0.000357 - 0.153290 * (sdss_i - sdss_z) - 0.009277 * (sdss_i - sdss_z)**2

    return { 'hsc_g': hsc_g, 'hsc_r': hsc_r, 'hsc_i': hsc_i }

def dm_to_kpc(dm):
    return 1e-3 * 10**((dm + 5.0) / 5.0)

def kpc_to_dm(kpc):
    return 5 * np.log(kpc * 1e3) / np.log(10) - 5

def arcmin_to_kpc(s, r):
    return 