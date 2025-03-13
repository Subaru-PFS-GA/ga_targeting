import logging
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time

from ..util.args import *
from ..util.coords import *
from ..util import safe_deep_copy, ReadOnlyDict
from ..photometry import Photometry, Color, Magnitude
from ..diagram import ColorAxis, MagnitudeAxis, DiagramValueProvider

from ..setup_logger import logger

class Catalog(DiagramValueProvider):
    def __init__(self, name=None, frame='icrs', equinox=None, epoch=None, orig=None):
        if not isinstance(orig, Catalog):
            self.__name = name
            self.__photometry = {}
            self.__frame = frame
            self.__equinox = equinox
            self.__epoch = epoch
        else:
            self.__name = name or orig.__name
            self.__photometry = orig.__photometry
            self.__frame = frame or orig.__frame
            self.__equinox = equinox or orig.__equinox
            self.__epoch = epoch or orig.__epoch

    def __len__(self):
        # TODO: return number of simulated stars per population
        raise NotImplementedError()

    def __get_name(self) -> str:
        return self.__name
    
    def __set_name(self, value: str):
        self.__name = value

    name = property(__get_name, __set_name)

    def __get_photometry(self):
        return ReadOnlyDict(self.__photometry)
    
    def _set_photometry(self, value):
        self.__photometry = value
    
    photometry = property(__get_photometry)

    def __get_frame(self):
        return self.__frame

    def __set_frame(self, value):
        self.__frame = value

    frame = property(__get_frame, __set_frame)

    def __get_equinox(self):
        return self.__equinox
    
    def __set_equinox(self, value):
        self.__equinox = value

    equinox = property(__get_equinox, __set_equinox)

    def __get_epoch(self):
        return self.__epoch
    
    def __set_epoch(self, value):
        self.__epoch = value

    epoch = property(__get_epoch, __set_epoch)

    def __get_observed(self):
        raise NotImplementedError()

    observed = property(__get_observed)

    def append_photometry(self, photometry):
        self.__photometry[photometry.name] = photometry

    def get_magnitude(self, magnitude: Magnitude, observed=False, dered=None, mask=None):
        raise NotImplementedError()

    def get_color(self, color: Color, observed=False, dered=None, mask=None):
        m1, s_m1 = self.get_magnitude(color.magnitudes[0], observed=observed, dered=dered, mask=mask)
        m2, s_m2 = self.get_magnitude(color.magnitudes[1], observed=observed, dered=dered, mask=mask)

        # TODO: correlated error?
        if s_m1 is not None and s_m2 is not None:
            s_c = np.sqrt(s_m1**2 + s_m2**2)
        else:
            s_c = None

        return m1 - m2, s_c

    def get_id(self, id=None, mask=None):
        mask = np.s_[()] if mask is None else mask

        id = np.array(self.data[id or 'objid'])[mask]
        return id
    
    def has_coords(self, ra_column=None, dec_column=None):
        ra_column = ra_column if ra_column is not None else 'RA'
        dec_column = dec_column if dec_column is not None else 'Dec'

        return ra_column in self.data and dec_column in self.data

    def get_coords(self, ra_column=None, dec_column=None, mask=None, ctype=None):
        """
        Returns the coordinates in any of the well-known array formats from the catalog.
        If `ra_column` or `dec_column` are specified, it assumes that those columns contain
        the coordinates.
        """

        ra_column = ra_column if ra_column is not None else 'RA'
        dec_column = dec_column if dec_column is not None else 'Dec'
        
        mask = np.s_[()] if mask is None else mask

        ra = np.array(self.data[ra_column])[mask]
        dec = np.array(self.data[dec_column])[mask]

        if ctype is not None:
            _, coords = normalize_coords(ra, dec)
            return denormalize_coords(ctype, coords)
        else:
            return ra, dec

    def get_skycoords(self, mask=None, frame=None, equinox=None, include_pm=False, **kwargs):
        """
        Returns the coordinates as astropy SkyCoords. The coordinate frame and equinox are
        taken from the catalog by default but can be overriden by the function parameters.
        """

        # TODO: what if we have coordinate errors?

        mask = np.s_[:] if mask is None else mask
        frame = normalize_frame(frame or self.__frame, equinox or self.__equinox)
        ra, dec = self.get_coords(mask=mask)

        if not include_pm:
            coords = SkyCoord(ra, dec, unit=u.degree, frame=frame, **kwargs)
        else:
            pmra = np.array(self.data['pmra'][mask])
            pmdec = np.array(self.data['pmdec'][mask])

            # If the epoch is specified on class level, it overrides the data
            if self.__epoch is not None:
                epoch = np.full_like(ra, normalize_epoch(self.__epoch).jyear)
            else:
                epoch = np.array(self.data['epoch'][mask].apply(normalize_epoch).apply(lambda e: e.jyear))

            if 'parallax' in self.data:
                parallax = np.array(self.data['parallax'][mask])
            else:
                parallax = np.zeros_like(ra)
            
            # Assume RV = 0 when not available
            if 'rv' in self.data:
                rv = np.array(self.data['rv'][mask])
            else:
                rv = np.zeros_like(ra)

            # Sanitize the parallax, apply_space_motion requires a non-zero distance
            parallax[np.isnan(parallax) | (parallax < 1e-7)] = 1e-7

            # Sanitize the peculiar motion components
            pmra[np.isnan(pmra)] = 0.0
            pmdec[np.isnan(pmdec)] = 0.0
            rv[np.isnan(rv)] = 0.0
            
            coords = SkyCoord(
                ra = ra * u.deg,
                dec = dec * u.deg,
                distance = Distance(parallax=parallax * u.mas, allow_negative=False),
                pm_ra_cosdec = pmra * u.mas / u.yr,
                pm_dec = pmdec * u.mas / u.yr,
                radial_velocity = rv * u.km / u.s,
                obstime = Time(epoch, format='jyear'),
                frame = self.frame,
                equinox = self.equinox)
        
        return coords

    def transform_coords(self, target_frame='icrs', target_equinox=None, target_epoch=None, apply_motion=True, mask=None):
        """
        Converts the coordinates of the catalog to the specified frame. The input frame
        is taken from the catalog.
        """
        
        mask = mask if mask is not None else np.s_[:]

        # Compare the source and target frames to see if any conversion is necessary
        source_frame = normalize_frame(self.__frame, self.__equinox)
        target_frame = normalize_frame(target_frame, equinox=target_equinox)

        source_equinox = source_frame.equinox.value if hasattr(source_frame, 'equinox') else None
        target_equinox = target_frame.equinox.value if hasattr(target_frame, 'equinox') else None

        need_frame_conv = (source_frame.name.lower() != target_frame.name.lower() or
                           source_equinox != target_equinox)

        # Compare the source and target epochs to see if any spatial motion propagation
        # is necessary. Here we assume, that `self.__epoch` always overrides the epoch values in
        # the DataFrame, unless None. This is fine, because when `self.__epoch` is None, the values
        # from the DataFrame must be taken and we can't assume they're all the same
        source_epoch = normalize_epoch(self.__epoch)
        target_epoch = normalize_epoch(target_epoch)

        need_space_motion = (apply_motion and target_epoch is not None and 
                             (source_epoch is None or
                              source_epoch is not None and source_epoch.jyear != target_epoch.jyear))

        logger.info(f'Catalog: `{self.name}`, frame {source_frame.name} ({source_equinox}), epoch {source_epoch}, '
                    f'target frame {target_frame.name} ({target_equinox}), epoch {target_epoch}')

        if not need_frame_conv and not need_space_motion:
            logger.info(f'Source and target frames, coordinate epochs are the same for catalog `{self.name}`, '
                        f'no coordinate conversion necessary.')
            return

        # If proper motions are available, include them in the conversions
        pm_avail = 'pmra' in self.data and 'pmdec' in self.data and 'epoch' in self.data
        source_coords = self.get_skycoords(mask=mask, include_pm=pm_avail)

        # Check if we can and need to apply spatial motion
        # If so, this is done in the original coordinate frame of the catalog

        if not pm_avail and apply_motion:
            logger.warning(f'Proper motion not available in catalog `{self.name}` to apply space motion.')
        elif target_epoch is None and apply_motion:
            logger.warning(f'Target epoch for catalog `{self.name}` is not defined. Cannot apply space motion.')
        elif pm_avail and need_space_motion:
            logger.info(f'Applying space motion for catalog `{self.name}` to epoch J{target_epoch.jyear:0.2f}.')

            # Here we silently assume that the epoch is always given in Julian years as opposed to
            # decimal Gregorian years. GAIA uses Julian years for the coordinate epochs
            source_coords = source_coords.apply_space_motion(target_epoch)
            self.__epoch = target_epoch.jyear

        if not need_frame_conv:
            logger.info(f'Source and target frames for catalog `{self.name}` are the same. No conversion necessary. '
                        f'Keeping {target_frame.name} with equinox {target_equinox}')
            target_coords = source_coords
        else:
            logger.info(f'Source and target frames for catalog `{self.name}` are not the same. '
                        f'Converting from {source_frame.name} with equinox {source_equinox} '
                        f'to {target_frame.name} with equinox {target_equinox}.')

            target_coords = source_coords.transform_to(frame=target_frame)

        # Save the new coordinates to the DataFrame
        self.data.loc[mask, 'RA'] = target_coords.ra.degree
        self.data.loc[mask, 'Dec'] = target_coords.dec.degree

        if pm_avail:
            self.data.loc[mask, 'parallax'] = target_coords.distance.parallax.mas
            self.data.loc[mask, 'pmra'] = target_coords.pm_ra_cosdec.value
            self.data.loc[mask, 'pmdec'] = target_coords.pm_dec.value
            self.data.loc[mask, 'rv'] = target_coords.radial_velocity.value

            if target_epoch is not None and 'epoch' in self.data:
                self.data.loc[mask, 'epoch'] = target_epoch.jyear
    
    def cone_search(self, pos, rad):
        pos = normalize_pos(pos)
        rad = normalize_angle(rad, u.arcmin)

        coords = self.get_skycoords()
        d = pos.separation(coords)
        mask = (d.arcmin <= rad.arcmin)

        return mask
    
    def random_sample(self, p):
        n = len(self)
        c = np.random.choice(np.arange(n), int(n * p), replace=False)
        mask = np.full((n,), False)
        mask[c] = True

        return mask
    
    def cross_match(self, other, max_separation=1, mask=None, mask_other=None):
        """
        Cross match two catalogs.
        
        Use built-in astropy functionality to perform the matching.
        """

        max_separation = normalize_angle(max_separation, u.arcsec, allow_none=True)

        c1 = self.get_skycoords(mask=mask)
        c2 = other.get_skycoords(mask=mask_other)

        idx, separation, _ = c1.match_to_catalog_sky(c2)
        
        # Apply cut on separation
        if max_separation is not None:
            idx[separation > max_separation] = -1

        return idx, separation