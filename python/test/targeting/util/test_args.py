import os
from datetime import datetime
import numpy as np
from astropy import units as u
from astropy.units import Quantity
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord
from astropy.units.core import UnitConversionError

from test_base import TestBase
from pfs.ga.targeting.util.args import *

class UtilArgsTest(TestBase):
    def test_normalize_time(self):
        self.assertIsNone(normalize_time(None))
        self.assertIsNone(normalize_time(None, allow_none=True))
        self.assertRaises(TypeError, normalize_time, None, allow_none=False)

        self.assertIsInstance(normalize_time(Time(datetime.utcnow())), Time)
        self.assertIsInstance(normalize_time(datetime.utcnow()), Time)

    def test_normalize_angle(self):
        self.assertIsNone(normalize_angle(None))
        self.assertIsNone(normalize_angle(None, allow_none=True))
        self.assertRaises(TypeError, normalize_angle, None, allow_none=False)

        self.assertIsInstance(normalize_angle(Angle(20, u.degree)), Angle)
        self.assertIsInstance(normalize_angle(20), Angle)

    def test_normalize_epoch(self):
        self.assertIsNone(normalize_epoch(None))
        self.assertIsNone(normalize_epoch(None, allow_none=True))
        self.assertRaises(TypeError, normalize_epoch, None, allow_none=False)

        self.assertEqual(2000.0, normalize_epoch('J2000').jyear)
        self.assertEqual(2015.15, normalize_epoch('j2015.15').jyear)
        self.assertEqual(2015.0, normalize_epoch(2015.0).jyear)
        self.assertEqual(2015.0, normalize_epoch(2015).jyear)

    def test_normalize_pos(self):
        self.assertIsNone(normalize_pos(None))
        self.assertIsNone(normalize_pos(None, allow_none=True))
        self.assertRaises(TypeError, normalize_pos, None, allow_none=False)

        # Scalar

        self.assertIsInstance(normalize_pos(SkyCoord(10 * u.degree, 10 * u.degree)), SkyCoord)
        self.assertIsInstance(normalize_pos(10, 10), SkyCoord)
        self.assertIsInstance(normalize_pos([10, 10]), SkyCoord)
        self.assertIsInstance(normalize_pos([Angle(10 * u.degree), Angle(10 * u.degree)]), SkyCoord)
        self.assertIsInstance(normalize_pos(Angle([10, 10] * u.degree)), SkyCoord)
        self.assertIsInstance(normalize_pos([ '15h 09m 08.5s', '+67d 13m 21s' ]), SkyCoord)

        # Array
        self.assertIsInstance(normalize_pos(SkyCoord([10, 20] * u.degree, [10, 20] * u.degree)), SkyCoord)
        self.assertIsInstance(normalize_pos([10], [10]), SkyCoord)
        self.assertIsInstance(normalize_pos([10, 20], [10, 20]), SkyCoord)
        self.assertIsInstance(normalize_pos([[10, 20], [10, 20]]), SkyCoord)
        self.assertIsInstance(normalize_pos(np.array([[10, 20], [10, 20], [10, 20]])), SkyCoord)

    def test_normalize_pm(self):
        self.assertEqual((None, None, None, None), normalize_pm(None))
        self.assertEqual((None, None, None, None), normalize_pm(None, allow_none=True))
        self.assertRaises(TypeError, normalize_pm, None, allow_none=False)

        self.assertEqual((Quantity(0.1 * u.mas / u.yr), Quantity(0.1 * u.mas / u.yr), None, None), normalize_pm([0.1, 0.1]))
        self.assertEqual((Quantity(0.1 * u.mas / u.yr), Quantity(0.1 * u.mas / u.yr), None, None), normalize_pm([Quantity(0.1 * u.mas / u.yr), Quantity(0.1 * u.mas / u.yr)]))

        self.assertEqual((Quantity(0.1 * u.mas / u.yr),
                          Quantity(0.1 * u.mas / u.yr),
                          Quantity(0.01 * u.mas / u.yr),
                          Quantity(0.01 * u.mas / u.yr)),
                          normalize_pm([0.1, 0.1], [0.01, 0.01]))
    
    def test_normalize_RV(self):
        self.assertEqual((None, None), normalize_RV(None))
        self.assertEqual((None, None), normalize_RV(None, allow_none=True))
        self.assertRaises(TypeError, normalize_RV, None, allow_none=False)

        self.assertEqual((Quantity(100 * u.km / u.s), None), normalize_RV(100))
        self.assertEqual((Quantity(100 * u.km / u.s), None), normalize_RV(Quantity(100 * u.km / u.s)))

        self.assertEqual((Quantity(100 * u.km / u.s),
                          Quantity(1 * u.km / u.s)),
                          normalize_RV(100, 1.0))
    
    def test_normalize_DM(self):
        self.assertEqual((None, None), normalize_DM(None))
        self.assertEqual((None, None), normalize_DM(None, allow_none=True))
        self.assertRaises(TypeError, normalize_DM, None, allow_none=False)

        self.assertEqual((Quantity(10 * u.mag), None), normalize_DM(10))
        self.assertEqual((Quantity(10 * u.mag), None), normalize_DM(Quantity(10 * u.mag)))

        self.assertEqual((Quantity(10 * u.mag),
                          Quantity(0.1 * u.mag)),
                          normalize_DM(10, 0.1))
        
    def test_normalize_dist(self):
        self.assertEqual((None, None), normalize_dist(None))
        self.assertEqual((None, None), normalize_dist(None, allow_none=True))
        self.assertRaises(TypeError, normalize_dist, None, allow_none=False)

        self.assertEqual((Quantity(10 * u.kpc), None), normalize_dist(10))
        self.assertEqual((Quantity(10 * u.kpc), None), normalize_dist(Quantity(10 * u.kpc)))

        self.assertEqual((Quantity(10 * u.kpc),
                          Quantity(0.1 * u.kpc)),
                          normalize_dist(10, 0.1))
        
    def test_normalize_parallax(self):
        self.assertEqual((None, None), normalize_parallax(None))
        self.assertEqual((None, None), normalize_parallax(None, allow_none=True))
        self.assertRaises(TypeError, normalize_parallax, None, allow_none=False)

        self.assertEqual((Angle(10 * u.arcsec), None), normalize_parallax(10))
        self.assertEqual((Angle(10 * u.arcsec), None), normalize_parallax(Angle(10 * u.arcsec)))

        self.assertEqual((Angle(10 * u.arcsec),
                          Angle(0.1 * u.arcsec)),
                          normalize_parallax(10, 0.1))
        
    def test_normalize_exp_time(self):
        self.assertEqual(None, normalize_exp_time(None))
        self.assertEqual(None, normalize_exp_time(None, allow_none=True))
        self.assertRaises(TypeError, normalize_exp_time, None, allow_none=False)

        self.assertEqual(Quantity(10 * u.s), normalize_exp_time(10))
        self.assertEqual(Quantity(10 * u.s), normalize_exp_time(Quantity(10 * u.s)))
        self.assertEqual(Quantity(60 * u.s), normalize_exp_time(Quantity(1 * u.min)))
        
        self.assertRaises(UnitConversionError, normalize_exp_time, Quantity(1 * u.km))

    def test_normalize_equinox(self):
        eq = normalize_equinox(None)
        self.assertIsNone(eq)

        eq = normalize_equinox(None, allow_none=True)
        self.assertIsNone(eq)

        self.assertRaises(TypeError, lambda: normalize_equinox(None, allow_none=False))

        eq = normalize_equinox('J2000.0')
        self.assertEqual(2000.0, eq.jyear)

        eq = normalize_equinox('2000.0')
        self.assertEqual(2000.0, eq.jyear)

        eq = normalize_equinox(2000.0)
        self.assertEqual(2000.0, eq.jyear)

        eq = normalize_equinox(Time(2000, format='jyear'))
        self.assertEqual(2000.0, eq.jyear)

    def test_normalize_frame(self):
        frame = SkyCoord(0 * u.deg, 0 * u.deg)
        f = normalize_frame(frame)
        self.assertIsInstance(f, BaseRADecFrame)
        self.assertEqual('icrs', f.name.lower())
        self.assertRaises(AttributeError, lambda: f.equinox)

        frame = SkyCoord(0 * u.deg, 0 * u.deg, frame='fk5')
        f = normalize_frame(frame)
        self.assertIsInstance(f, BaseRADecFrame)
        self.assertEqual('fk5', f.name.lower())
        self.assertEqual(2000.0, f.equinox.jyear)

        frame = SkyCoord(0 * u.deg, 0 * u.deg, frame='fk5')
        f = normalize_frame(frame, equinox=2015.0)
        self.assertIsInstance(f, BaseRADecFrame)
        self.assertEqual('fk5', f.name.lower())
        self.assertEqual(2015.0, f.equinox.jyear)
        
        f = normalize_frame(f)
        self.assertIsInstance(f, BaseRADecFrame)
        self.assertEqual('fk5', f.name.lower())
        self.assertEqual(2015.0, f.equinox.jyear)

        frame, equinox = 'fk5', 2000.0
        f = normalize_frame(frame, equinox=equinox)
        self.assertIsInstance(f, BaseRADecFrame)
        self.assertEqual('fk5', f.name.lower())
        self.assertEqual(2000.0, f.equinox.jyear)

    def test_normalize_skycoord(self):
        pos, pm_err, RV_err, obs_time = normalize_skycoord(12, 14, pm=[2, 2], pm_err=[0.2, 0.2], RV=125, RV_err=3)
        pos, pm_err, RV_err, obs_time = normalize_skycoord([12, 14], pm=[2, 2], pm_err=[0.2, 0.2], RV=125, RV_err=3)

    def test_normalize_distance(self):
        distance, dist_err, DM_err, parallax_err = normalize_distance(dist=10, dist_err=0.1)
        self.assertEqual(distance.value, 10)

        distance, dist_err, DM_err, parallax_err = normalize_distance(DM=10, DM_err=0.1)
        self.assertEqual(distance.distmod.value, 10)

        distance, dist_err, DM_err, parallax_err = normalize_distance(parallax=10, parallax_err=0.1)
        self.assertAlmostEqual(distance.parallax.arcsec, 10)
        self.assertAlmostEqual(parallax_err.arcsec, 0.1)