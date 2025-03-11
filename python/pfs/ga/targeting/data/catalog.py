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
    
    def has_coords(self, ra=None, dec=None):
        ra = ra if ra is not None else 'RA'
        dec = dec if dec is not None else 'Dec'

        return ra in self.data and dec in self.data

    def get_coords(self, ra_column=None, dec=None, mask=None, ctype=None):
        """
        Returns the coordinates in any of the well-known array formats from the catalog.
        If `ra_column` or `dec_column` are specified, it assumes that those columns contain
        the coordinates.
        """

        ra_column = ra_column if ra_column is not None else 'RA'
        dec = dec if dec is not None else 'Dec'
        
        mask = np.s_[()] if mask is None else mask

        ra_column = np.array(self.data[ra_column])[mask]
        dec = np.array(self.data[dec])[mask]

        if ctype is not None:
            _, coords = normalize_coords(ra_column, dec)
            return denormalize_coords(ctype, coords)
        else:
            return ra_column, dec

    def get_skycoords(self, mask=None, frame=None, equinox=None, **kwargs):
        """
        Returns the coordinates as astropy SkyCoords. The coordinate frame and equinox are
        taken from the catalog by default but can be overriden by the function parameters.
        """

        # TODO: what if we have coordinate errors?

        # TODO: include parallax, etc, if available

        frame = frame or self.__frame
        equinox = equinox or self.__equinox

        ra, dec = self.get_coords(mask=mask)
        coords = SkyCoord(ra, dec, unit=u.degree, frame=frame, equinox=equinox, **kwargs)
        return coords

    def transform_coords(self, target_frame='icrs', target_equinox=None, target_epoch=None):
        """
        Converts the coordinates of the catalog to the specified frame. The input frame
        is taken from the catalog.
        """

        # Create a target frame here, with an optional equinox, because `transform_to` doesn't allow for
        # specifying the equinox for FK5 etc.
        if target_equinox is None:
            logger.info(f'Target frame is {target_frame}.')
            target_frame = SkyCoord(0 * u.deg, 0 * u.deg, frame=target_frame).frame
        else:
            logger.info(f'Target frame is {target_frame}, target equinox is {target_equinox}.')
            target_frame = SkyCoord(0 * u.deg, 0 * u.deg, frame=target_frame, equinox=Time(target_equinox, format='decimalyear')).frame

        if 'pmra' in self.data and 'pmdec' in self.data:
            logger.info(f'Transforming coordinates and proper motions.')
            self.__transform_coords_pm(target_frame=target_frame, target_epoch=target_epoch)
        else:
            logger.info(f'Transforming coordinates only, without proper motion.')
            self.__transform_coords(target_frame=target_frame)

    def __transform_coords(self, target_frame='icrs'):
        ra = np.array(self.data['RA'])
        dec = np.array(self.data['Dec'])

        # Coordinates in the source catalog
        source_coords = SkyCoord(
            ra = ra * u.deg,
            dec = dec * u.deg,
            frame = self.frame,
            equinox = self.equinox)

        # Convert coordinates to the default target frame
        target_coords = source_coords.transform_to(frame=target_frame)
        
        self.data['RA'] = target_coords.ra.degree
        self.data['Dec'] = target_coords.dec.degree

    def __transform_coords_pm(self, target_frame='icrs', target_epoch=None):
        """
        Transform the coordinates into a new frame and propagate them to a new epoch if
        proper motions are available.
        """

        # Astropy only takes a single epoch when applying proper motions, so find all unique epochs
        # and apply proper motion in batches
        unique_epochs = self.data['epoch'].unique()

        logger.info(f'{unique_epochs.size} unique epochs found in the catalog.')

        for epoch in unique_epochs:
            # Filter down the dataset to this epoch only
            mask = (self.data['epoch'] == epoch)
            epoch = normalize_epoch(epoch)

            # Skip coordinate conversion when frame and epoch are the same
            # if (epoch == target_epoch
            #     and catalog.frame.lower() == target_frame.lower()
            #     and catalog.equinox == target_equinox):

            #     continue

            ra = np.array(self.data['RA'][mask])
            dec = np.array(self.data['Dec'][mask])
            pmra = np.array(self.data['pmra'][mask])
            pmdec = np.array(self.data['pmdec'][mask])

            if 'parallax' in self.data:
                parallax = np.array(self.data['parallax'][mask])
            else:
                parallax = np.zeros_like(ra)
            
            # Assume RV = 0 when not available
            if 'rv' in self.data:
                rv = np.array(self.data['rv'][mask])
            else:
                rv = np.zeros_like(ra)

            # Sanitize the peculiar motion components
            parallax[np.isnan(parallax) | (parallax < 1e-7)] = 1e-7
            pmra[np.isnan(pmra)] = 0.0
            pmdec[np.isnan(pmdec)] = 0.0
            rv[np.isnan(rv)] = 0.0

            # Coordinates in the source catalog
            source_coords = SkyCoord(
                ra = ra * u.deg,
                dec = dec * u.deg,
                distance = Distance(parallax=parallax * u.mas, allow_negative=False),
                pm_ra_cosdec = pmra * u.mas / u.yr,
                pm_dec = pmdec * u.mas / u.yr,
                radial_velocity = rv * u.km / u.s,
                obstime = Time(epoch, format='decimalyear'),
                frame = self.frame,
                equinox = self.equinox)

            # Propagate coordinates to the default epoch
            logger.info(f'Propagating coordinates from epoch {epoch} to target epoch {target_epoch}.')
            target_coords = source_coords.apply_space_motion(Time(target_epoch, format='decimalyear'))

            # Convert coordinates to the default target frame
            logger.info(f'Transforming coordinates to target frame.')
            target_coords = target_coords.transform_to(frame=target_frame)

            self.data.loc[mask, 'RA'] = target_coords.ra.degree
            self.data.loc[mask, 'Dec'] = target_coords.dec.degree
            self.data.loc[mask, 'parallax'] = target_coords.distance.parallax.mas
            self.data.loc[mask, 'pmra'] = target_coords.pm_ra_cosdec.value
            self.data.loc[mask, 'pmdec'] = target_coords.pm_dec.value
            self.data.loc[mask, 'rv'] = target_coords.radial_velocity.value
            if self.data['epoch'].dtype == 'string' or self.data['epoch'].dtype == str:
                self.data.loc[mask, 'epoch'] = f'J{target_epoch}'
            else:
                self.data.loc[mask, 'epoch'] = target_epoch
    
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