from collections.abc import Iterable
from numbers import Number
from datetime import datetime
import numpy as np
from astropy import units as u
from astropy.units import Quantity
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, Distance, BaseRADecFrame

def __unwrap_iterable(v):
    if v is None:
        return None
    if isinstance(v, (Angle, SkyCoord)):            # Need to test it because everything is iterable in astropy
        return v
    elif isinstance(v, Iterable) and len(v) == 1:
        return v[0]
    else:
        return v

def __is_n_of_type(v, type, n):
    if isinstance(v, Iterable):
        if len(v) == n:
            for i in range(n):
                if not isinstance(v[i], type):
                    return False
            return True
        else:
            return False
    else:
        return False
    
def __is_two_of_type(v, type):
    return __is_n_of_type(v, type, 2)

def normalize_array(a, allow_none=True):
    if allow_none and a is None:
        return None
    elif a is None:
        raise ValueError('Argument cannot be None.')

    if isinstance(a, np.ndarray):
        return a

    if not isinstance(a, Iterable):
        a = [a]

    return np.array(a)

def normalize_time(time, scale='utc', allow_none=True):
    if allow_none and time is None:
        pass
    elif time is None:
        raise TypeError('Argument `time` cannot be None.')
    elif isinstance(time, Time):
        pass
    elif isinstance(time, datetime):
        time = Time(time, format='datetime', scale=scale)
    elif isinstance(time, str):
        # Assume ISO time
        time = Time(time, scale=scale)
    else:
        raise NotImplementedError()
    
    return time

def normalize_angle(angle, unit=None, allow_none=True):
    # TODO: Add support for arrays

    if allow_none and angle is None:
        return None
    elif angle is None:
        raise TypeError('Argument `angle` cannot be None.')
    elif isinstance(angle, Angle):
        pass
    elif isinstance(angle, Quantity):
        angle = Angle(angle)
    elif isinstance(angle, str):
        angle = Angle(angle)
    elif isinstance(angle, Number):
        unit = unit if unit is not None else u.degree
        angle = Angle(np.float64(angle), unit=unit)
    else:
        raise NotImplementedError()
    
    return angle

def normalize_epoch(epoch, allow_none=True):
    if allow_none and epoch is None:
        return None
    elif epoch is None:
        raise TypeError('Argument `epoch` cannot be None.')
    elif isinstance(epoch, Time):
        pass
    elif isinstance(epoch, str):
        if epoch.upper().startswith('J'):
            epoch = Time(float(epoch[1:]), format='jyear')
        else:
            epoch = Time(float(epoch), format='jyear')
    elif isinstance(epoch, Number):
        epoch = Time(epoch, format='jyear')
    else:
        raise NotImplementedError()

    return epoch

def normalize_pos(*pos, unit=u.degree, frame='icrs', equinox=None, allow_none=True):
    # Bring a pair of coordinates into a normalized format

    pos = __unwrap_iterable(pos)

    if allow_none and pos is None:
        pass
    elif pos is None:
        raise TypeError('Argument `pos` cannot be None.')
    elif isinstance(pos, SkyCoord) and \
         pos.frame.name == frame and \
         (frame == 'icrs' or (pos.equinox is None or pos.equinox == equinox)):
        
        pass
    elif isinstance(pos, SkyCoord):
        # TODO: add frame conversion
        raise NotImplementedError()
    elif isinstance(pos, np.ndarray) and pos.shape[-1] == 2:
        pos = SkyCoord(np.float64(pos[..., 0]) * unit, np.float64(pos[..., 1]) * unit, frame=frame, equinox=equinox)
    elif isinstance(pos, Iterable):
        if __is_two_of_type(pos, Number):
            pos = SkyCoord(np.float64(pos[0]) * unit, np.float64(pos[1]) * unit, frame=frame, equinox=equinox)
        elif __is_two_of_type(pos, str):
            pos = SkyCoord(Angle(pos[0]), Angle(pos[1]), unit=u.deg, frame=frame, equinox=equinox)
        elif __is_two_of_type(pos, Angle):
            pos = SkyCoord(pos[0], pos[1], frame=frame, equinox=equinox)
        elif __is_two_of_type(pos, Iterable):
            pos = SkyCoord(np.float64(pos[0]) * unit, np.float64(pos[1]) * unit, frame=frame, equinox=equinox)
        else:
            raise NotImplementedError()
    else:
        raise NotImplementedError()
    
    return pos

def normalize_pm(pm, pm_err=None, allow_none=True):
    pm = __unwrap_iterable(pm)
    pm_err = __unwrap_iterable(pm_err)

    if allow_none and pm is None:
        pm_ra, pm_dec = None, None
    elif pm is None:
        raise TypeError('Argument `pm` cannot be None.')
    elif __is_two_of_type(pm, Quantity) and pm[0].unit == u.mas / u.yr and pm[1].unit == u.mas / u.yr:
        [ pm_ra, pm_dec ] = pm
    elif __is_two_of_type(pm, Quantity):
        raise NotImplementedError()
    elif __is_two_of_type(pm, Number):
        [ pm_ra, pm_dec ] = [ Quantity(np.float64(pm[0]), u.mas / u.yr), Quantity(np.float64(pm[1]), u.mas / u.yr) ]
    else:
        raise NotImplementedError()

    if pm_err is None:
        pm_ra_err, pm_dec_err = None, None
    elif __is_two_of_type(pm, Quantity) and pm_err.unit == u.mas / u.yr:
        [ pm_ra_err, pm_dec_err ] = pm_err
    elif __is_two_of_type(pm_err, Number):
        [ pm_ra_err, pm_dec_err ] = [ Quantity(np.float64(pm_err[0]), u.mas / u.yr), Quantity(np.float64(pm_err[1]), u.mas / u.yr) ]
    else:
        raise NotImplementedError()
    
    return pm_ra, pm_dec, pm_ra_err, pm_dec_err

def __normalize_scalar_quantity(name, q, q_err, unit, allow_none=True, cls=Quantity):
    if allow_none and q is None:
        pass
    elif q is None:
        raise TypeError(f'Argument {name}` cannot be None.')
    elif isinstance(q, cls) and q.unit == unit:
        pass
    elif isinstance(q, cls):
        # Try to convert to the desired unit
        q = q.to(unit)
    elif isinstance(q, Number):
        q = cls(np.float64(q), unit=unit)
    else:
        raise NotImplementedError()
    
    if q_err is None:
        pass
    elif isinstance(q_err, cls) and q_err.unit == unit:
        pass
    elif isinstance(q_err, cls):
        # Try to convert to the desired unit
        q = q.to(unit)
    elif isinstance(q_err, Number):
        q_err = cls(np.float64(q_err), unit)
    else:
        raise NotImplementedError()

    return q, q_err

def normalize_RV(RV, RV_err=None, allow_none=True):
    return __normalize_scalar_quantity('RV', RV, RV_err, u.km / u.s, allow_none=allow_none)

def normalize_DM(DM, DM_err=None, allow_none=True):
    return __normalize_scalar_quantity('DM', DM, DM_err, u.mag, allow_none=allow_none)

def normalize_dist(dist, dist_err=None, allow_none=True):
    return __normalize_scalar_quantity('dist', dist, dist_err, u.kpc, allow_none=allow_none)

def normalize_parallax(parallax, parallax_err=None, allow_none=True):
    return __normalize_scalar_quantity('parallax', parallax, parallax_err, u.arcsec, allow_none=allow_none, cls=Angle)

def normalize_exp_time(exp_time, allow_none=True):
    exp_time, _ = __normalize_scalar_quantity('exp_time', exp_time, None, u.s, allow_none=allow_none)
    return exp_time

def normalize_equinox(equinox, allow_none=True):
    if allow_none and equinox is None:
        return None
    elif equinox is None:
        raise TypeError(f'Argument {equinox}` cannot be None.')
    elif isinstance(equinox, Time):
        pass
    elif isinstance(equinox, str):
        if equinox.lower().startswith('j'):
            equinox = Time(equinox, format='jyear_str')
        else:
            equinox = Time(float(equinox), format='jyear')
    else:
        equinox = Time(equinox, format='jyear')

    return equinox

def normalize_frame(frame, equinox=None):
    """
    Normalize a frame definition.

    Parameters
    ----------
    frame : string or SkyCoord or Frame
        String abbreviation of the frame or a SkyCoord defined in the frame.
    equinox : Time or float
        Equinox of the frame, if applicable. Will be ignored if `frame` is SkyCoord.
    """

    equinox = normalize_equinox(equinox, allow_none=True)

    # Create a target frame through the SkyCoord initializer, with an optional equinox, because
    # `transform_to` doesn't allow for specifying the equinox for FK5 etc.

    if isinstance(frame, SkyCoord):
        if equinox is None:
            # Just return the frame information from the SkyCoord isntance
            frame = frame.frame
        else:
            # Create a new frame with the specified equinox
            frame = SkyCoord(0 * u.deg, 0 * u.deg, frame=type(frame.frame), equinox=equinox).frame
    elif isinstance(frame, BaseRADecFrame):
        if equinox is None:
            pass
        else:
            frame = SkyCoord(0 * u.deg, 0 * u.deg, frame=type(frame), equinox=equinox).frame
    else:
        frame = SkyCoord(0 * u.deg, 0 * u.deg, frame=frame, equinox=equinox).frame
        
    return frame

def normalize_skycoord(*pos, pm=None, pm_err=None, RV=None, RV_err=None, obs_time=None, frame='icrs', equinox=None, scale='utc'):
    # TODO: Add support for arrays

    pos = normalize_pos(*pos, frame=frame, equinox=equinox)
    pm_ra, pm_dec, pm_ra_err, pm_dec_err = normalize_pm(pm, pm_err)
    RV, RV_err = normalize_RV(RV, RV_err)
    obs_time = normalize_time(obs_time, scale=scale)

    if isinstance(pos, SkyCoord) and pm is None and RV is None and obs_time is None and \
       (frame is None or pos.frame.name == frame) and (equinox is None or pos.equinox == equinox):
        pass
    else:
        pos = SkyCoord(pos.ra, pos.dec,
                    pm_ra_cosdec=pm_ra, pm_dec=pm_dec,
                    radial_velocity=RV,
                    obstime=obs_time,
                    frame=frame, equinox=equinox)
    
    if pm_err is not None:
        pm_err = np.stack([ pm_ra_err, pm_dec_err ])

    return pos, pm_err, RV_err, obs_time

def normalize_distance(dist=None, dist_err=None, DM=None, DM_err=None, parallax=None, parallax_err=None, allow_none=True):
    dist, dist_err = normalize_dist(dist, dist_err)
    DM, DM_err = normalize_DM(DM, DM_err)
    parallax, parallax_err = normalize_parallax(parallax, parallax_err)

    if allow_none and dist is None and DM is None and parallax is None:
        distance = None
    elif (dist is not None or DM is not None) and parallax is None:
        distance = Distance(value=dist, distmod=DM)
    elif parallax is not None:
        distance = Distance(parallax=parallax)
    else:
        raise NotImplementedError()
    
    return distance, dist_err, DM_err, parallax_err