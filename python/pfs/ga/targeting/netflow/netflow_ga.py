#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 15:11:40 2023

@author: elenasaez, evankirby
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dill as pickle
import os
import h5py

from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u

import ets_fiber_assigner.netflow as nf
#from ics.cobraOps.TargetGroup import TargetGroup
#from ics.cobraOps.CollisionSimulator import CollisionSimulator
#from ics.cobraOps.CollisionSimulator2 import CollisionSimulator2
from ics.cobraOps.Bench import Bench
#from ics.cobraOps.cobraConstants import NULL_TARGET_POSITION, NULL_TARGET_ID
#from ics.cobraOps import plotUtils
from collections import defaultdict

from astropy.table import Table, vstack

import gurobipy as gbp

class GA_Netflow:

    def __init__(self, galaxy_name='Fornax', read_pickle=False, first_pointing_only=True):
        galaxy = self.Galaxy(galaxy_name)
        
        problemfile = galaxy_name.lower() + '_netflow_problem.mps'
        solutionfile = galaxy_name.lower() + 'netflow_problem.sol'
        picklefile1 = galaxy_name.lower() + '_netflow_galaxy.pkl'
        picklefile2 = galaxy_name.lower() + '_netflow_targets.pkl'
        picklefile3 = galaxy_name.lower() + '_netflow_fiber_assignments.pkl'
        
        if read_pickle and os.path.exists(picklefile1) and os.path.exists(picklefile2):
            #problem = gbp.read(solutionfile)
            #problem.optimize()
            with open(picklefile1, 'rb') as file:
                catalog = pickle.load(file)
            with open(picklefile2, 'rb') as file:
                targets = pickle.load(file)
        else:
            #Assign priorities to each star in the catalog file
            catalog, targets = self.assign_priorities(galaxy, use_prob=galaxy_name.lower() != 'm31')
            self.plot_netflow(catalog, plot_selected_only=False)

            with open(picklefile1, 'wb') as file:
                pickle.dump(catalog, file)
            with open(picklefile2, 'wb') as file:
                pickle.dump(targets, file)

        if read_pickle and os.path.exists(picklefile3):
            with open(picklefile3, 'rb') as file:
                res = pickle.load(file)
                fp_pos = pickle.load(file)
                telescopes = pickle.load(file)        
                
        if first_pointing_only:
            catalog, targets = self.run_netflow(galaxy, catalog, targets, ra=galaxy.ra0[0], dec=galaxy.dec0[0], pa=galaxy.pa[0], picklefile=picklefile3, read_pickle=read_pickle)
        else:
            for i, (r, d, p) in enumerate(zip(galaxy.ra0, galaxy.dec0, galaxy.pa)):
                catalog, targets = self.run_netflow(galaxy, catalog, targets, ra=r, dec=d, pa=p, picklefile=picklefile3, read_pickle=read_pickle)
                
        self.count_targets(catalog)
        
        #with open(picklefile1, 'wb') as file:
            #pickle.dump(problem, file)
        with open(picklefile1, 'wb') as file:
            pickle.dump(catalog, file)
        with open(picklefile2, 'wb') as file:
            pickle.dump(targets, file)
        
        self.plot_netflow(catalog, plot_selected_only=True)
        self.plot_hexagon(galaxy, catalog)
        
        
    class Galaxy:
         #Define the galaxy object including its ra, dec, file name, fluxstd file, and sky file
        def __init__(self, name):
            self.name = name
            if name.lower() != 'm31':
                self.photfile = 'hsc/'+name.lower()+'_tpall3e_g24.cat'
            self.probfile = 'pmap_'+name.lower()+'.h5'
            self.fluxstd = pd.read_feather('feather/fluxstd_'+name.lower()+'.feather')
            self.sky = pd.read_feather('feather/sky_'+name.lower()+'.feather')
            self.name = name
            if(name.lower() =='ursaminor'):
                """
                self.ra0 = [228.2,226.3,226.0,228.0]
                self.dec0 = [67.5,67.5,66.9,66.955]
                self.pa = [0,0,0,0]
                """
                self.ra = 227.29725 
                self.dec = 67.21436111
                self.ra0 = [self.ra]
                self.dec0 = [self.dec]
                self.pa = [0]
                self.pmra = -0.119     #Pace et al. (2022)
                self.pmdec = 0.072
                self.pmraerr = 0.005
                self.pmdecerr = 0.005
            elif(name.lower() == 'draco'):
                self.ra0 = [259.2,260.9,260.1,260.0]
                self.dec0 = [57.871,57.957,57.77,58.05]
                self.pa = [0,0,30,30]
                self.ra = 260.051666667
                self.dec = 57.9152777778
                self.pmra = 0.046
                self.pmdec = -0.188
                self.pmraerr = 0.006
                self.pmdecerr = 0.006
            elif (name.lower() == 'fornax'):
                self.ra0 = [39.5,40.4,40.3,39.5,40.9,40.5,39.4,39.1]
                self.dec0 = [-34.2,-34.1,-34.8,-34.9,-33.8,-35.1,-34.0,-35.2]
                self.pa = [30,30,30,30,30,30,30,30]
                self.ra = 39.9971
                self.dec = -34.4492
                self.pmra = 0.381
                self.pmdec = -0.358
                self.pmraerr = 0.001
                self.pmdecerr = 0.002
            elif(name.lower() == 'sculptor'):
                self.ra0 = [14.5,15.1,15.5,15.0,16.4,15.06,13.7,14.9]
                self.dec0 = [-33.7,-33.4,-33.7,-34.1,-33.9,-33.0,-33.55,-34.5]
                self.pa = [30,30,30,30,30,30,30,30]
                self.ra = 15.0392
                self.dec = -33.7089
                self.pmra = 0.101
                self.pmdec = -0.156
                self.pmraerr = 0.003
                self.pmdecerr = 0.002
            elif(name.lower() == 'bootes'):
                self.ra0 = [210.5,209.6,210.1,210.1]
                self.dec0 = [14.5,14.5,14.15,14.8]
                self.pa = [30,30,30,30]
                self.ra = 210.025
                self.dec = 14.5
                self.pmra = -0.387
                self.pmdec = -1.064
                self.pmraerr = 0.122
                self.pmdecerr = 0.098
            elif(name.lower() == 'ngc6822'):
                self.ra0 = [296.235]
                self.dec0 = [-14.789]
                self.pa = [30]
                self.ra = 296.234
                self.dec = -14.7976
            elif(name.lower() == 'm31'):
                ra_angle = Angle(['00:34:54', '00:28:56', '00:22:49', '00:41:43', '00:35:47', '00:29:42', '00:23:26', '00:42:45', '00:36:42', '00:37:40', '00:34:04', '00:28:13', '00:22:13', '00:50:18', '00:55:44', '01:01:02', '00:51:28', '00:56:59', '01:02:22', '00:52:41', '00:58:18', '01:03:46', '00:43:43', '00:49:12', '00:54:32', '00:37:13', '00:42:44', '00:48:06', '00:16:31.7', '00:10:04.7', '00:10:24.5', '00:16:04.2', '00:09:45.8', '00:06:10', '00:52:10', '00:48:10', '00:47:06', '00:47:05', '00:39:30', '00:37:40', '00:32:40', '00:31:10'], u.hour)
                dec_angle = Angle(['+42:22:57', '+43:11:01', '+43:57:54', '+42:54:25', '+43:43:47', '+44:32:01', '+45:19:01', '+44:15:02', '+45:04:35', '+46:25:20', '+41:02:05', '+41:50:00', '+42:36:45', '+40:07:25', '+39:15:27', '+38:22:34', '+41:27:45', '+40:35:34', '+39:42:28', '+42:48:01', '+41:55:37', '+41:02:17', '+39:37:51', '+38:47:03', '+37:55:17', '+39:06:52', '+38:17:16', '+37:26:39', '+44:43:30', '+45:27:47', '+46:49:07', '+43:22:15', '+44:06:28', '+47:35:07', '+44:25:07', '+43:30:07', '+42:16:07', '+40:40:07', '+41:50:07', '+40:05:04', '+40:05:04', '+39:05:04'], u.degree)
                self.ra0 = ra_angle.degree
                self.dec0 = dec_angle.degree
                self.pa = np.zeros(len(self.ra0))
                self.ra = 10.684708333
                self.dec = 41.26875
                self.photfile = 'hsc/M31Catalog_forPFS.csv'
            else:
                raise Exception("I don't know this galaxy.")
    
    
    def assign_probabilities(self, catalog, edges_x, edges_y, hist):
        gi = (catalog['gpsf']-catalog['ipsf'])
        gmag = catalog['gpsf']
        prob = np.zeros(len(catalog))
        inbounds = np.zeros(gi.shape,dtype=bool)
        inbounds[(gi >= np.min(edges_x)) & (gi < np.max(edges_x)) & (gmag >= np.min(edges_y)) & (gmag < np.max(edges_y))] = True
        prob[~inbounds] = 0.0
        
        x = np.digitize(gi[inbounds], edges_x, right=True) - 1
        y = np.digitize(gmag[inbounds], edges_y, right=False) - 1
        hist[0,x,y] /= np.sum(hist[0,x,y])
        hist[1,x,y] /= np.sum(hist[1,x,y])
        prob[inbounds] = hist[1,x,y] / (hist[0,x,y]+hist[1,x,y])
        
        prob[np.isnan(prob)] = 0.0
        
        return(prob)
    
    
    def get_gaia(self, galaxy, rad_arcmin=60):
        from astroquery.gaia import Gaia
        gbprp = ['G', 'BP', 'RP']
        Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaiadr3.gaia_source"
        Gaia.ROW_LIMIT = -1
        
        c_galaxy = SkyCoord(galaxy.ra*u.degree, galaxy.dec*u.degree)
        rastr = c_galaxy.ra.to_string(decimal=True, precision=8)
        decstr = c_galaxy.dec.to_string(decimal=True, precision=8)
        querystr = """
        SELECT source_id, ra, dec, parallax, parallax_error, pm, pmra, pmra_error, pmdec, pmdec_error, phot_g_mean_mag, phot_g_mean_flux_over_error, phot_bp_mean_mag, phot_bp_mean_flux_over_error, phot_rp_mean_mag, phot_rp_mean_flux_over_error
            FROM gaiadr3.gaia_source
            WHERE CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),
            CIRCLE('ICRS',{},{},{}))=1;""".format(rastr, decstr, str(rad_arcmin/60.))
        job = Gaia.launch_job_async(querystr, dump_to_file=False)
        gaia = job.get_results()
        
        gaia.rename_column('ra', 'ra_Gaia')
        gaia.rename_column('dec', 'dec_Gaia')
        for filter in gbprp:
            gaia.rename_column('phot_'+filter.lower()+'_mean_mag', filter)
            gaia[filter+'err'] = 2.5 / (np.log(10.) * gaia['phot_'+filter.lower()+'_mean_flux_over_error'])
            gaia.remove_column('phot_'+filter.lower()+'_mean_flux_over_error')
        gaia['coord'] = SkyCoord(gaia['ra_Gaia'], gaia['dec_Gaia'], unit=(u.degree, u.degree))
        
        # Include a column to show that this object was detected in Gaia.
        gaia['Gaia'] = True
        
        return gaia
    
    def merge_catalogs(self, cat1, cat2, cat_names=['1', '2'], match_radius=0.1, join_type='outer'):
        from astropy.table import join, join_skycoord, vstack
        
        """Make the intersection of the two catalogs."""
        join_func = join_skycoord(match_radius * u.arcsec)
        j = join_func(cat1['coord'], cat2['coord'])
        data = join(cat1, cat2, keys='coord', join_funcs={'coord': join_func}, 
                         table_names=cat_names)
        
        if (join_type == 'outer'):
            # Remove SkyCoords from table because outer join (i.e., vstack)
            # will not work with SkyCoords.
            for cat in [data, cat1, cat2]:
                removec = []
                for c in cat.columns:
                    if c.startswith('coord'):
                        removec.append(c)
                cat.remove_columns(removec)
            
            # Find the objects that are present in HSC but not Gaia.
            _, ind1, _ = np.intersect1d(j[0], np.setdiff1d(j[0], j[1]), 
                                        return_indices=True)
            data = vstack([data, cat1[ind1]])
            
            # Put SkyCoords back in the Table.
            cats = ['Gaia', 'HSC']
            for c in cats:
                if c in data.columns:
                    data['coord_'+c] = SkyCoord(data['ra_'+c], data['dec_'+c], unit=(u.degree, u.degree))
            
            # Define a primary coordinate, with primacy ordered by cats, defined above.
            for c in cats[::-1]:
                if c in data.columns:
                    coord = data['coord_'+c]
            for c in cats[::-1]:
                if c in data.columns:
                    coord[data[c]] = data['coord_'+c][data[c]]
            data['coord'] = coord
            
        return data
    
    def assign_priorities(self, galaxy, use_prob=True, showplot=True, write_catalog=True):
        if galaxy.name.lower() == 'm31':
            hsc = ascii.read(galaxy.photfile)
            hsc['id'] = np.arange(len(hsc), dtype=float)+1
        else:
            hsc = ascii.read(galaxy.photfile, 
                             names=['id', 'RA', 'Dec', 'X', 'Y', 'ipsf', 'gpsf', 'npsf', 
                                    'ipsferr', 'gpsferr', 'npsferr', 'cli', 'clg', 'cln', 
                                    'a_i', 'a_g', 'a_n'])
    
        hsc['coord'] = SkyCoord(hsc['RA'], hsc['Dec'], unit=(u.degree, u.degree))
        
        hsc['g0'] = hsc['gpsf'] - hsc['a_g']
        hsc['i0'] = hsc['ipsf'] - hsc['a_i']
        hsc['n0'] = hsc['npsf'] - hsc['a_n']
        
        hsc['id'] = hsc['id'].astype(float)
    
        hsc['selected'] = 0
        hsc['alreadyObserved'] = 0
        hsc['targetClass'] = 0
        hsc['exptime'] = 0
        
        hsc['priority'] = 0
        hsc['code'] = 0
        
        gaia = self.get_gaia(galaxy)
        catalog = self.merge_catalogs(hsc, gaia, cat_names=['HSC', 'Gaia'])
        
        g0 = catalog['g0']
        i0 = catalog['i0']
        gi0 = catalog['g0']-catalog['i0']
        
        targets = []
        
        if use_prob:
            def nb_cut(g, i, nb):
                return (0.12 < g - i) & (g - i < 0.5) | (g - nb > 0.1) & (g - i < 1.65) | (-0.25 * (g - i) + (g - nb) > -0.15)
            
            with h5py.File(galaxy.probfile, 'r') as f:
                edges_x = f['pmap']['edges_0'][:]
                edges_y = f['pmap']['edges_1'][:]
                hist = f['pmap']['hist'][:]
            
            catalog['prob'] = self.assign_probabilities(catalog, edges_x, edges_y, hist)
            catalog['prob'][~nb_cut(g0, i0, catalog['n0'])] = 0
            top_pri = np.maximum(np.floor((i0 - 16)/(21.5-16) * 8).astype(int)-7, -7) #top pri goes from 0-4 based on brightness 
            bot_pri = np.maximum(np.floor((i0 - 16)/(21.5-16) * 6).astype(int)+3, 3) #bot pri goes from 3-8 based on brightness 
            catalog['priority'] = np.minimum(np.maximum(bot_pri - np.rint(catalog['prob']*(bot_pri - top_pri)).astype(int), 0), 9)
                        
            w = (catalog['prob'] == 0.0)
            catalog['priority'][w] = 9
            catalog['code'][w] = 0
            
            w = (g0 > 19.6) & (g0 < 20.2) & (gi0 > -0.5) & (gi0 < 0.2) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 6
            catalog['code'][w] = 0
            
            w = (i0 <= 16) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 9
            catalog['code'][w] = 1
            
            w = (i0 >= 21.5) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 9
            catalog['code'][w] = 2
            
            w = (catalog['cli'] > 0.5) | (catalog['clg'] > 0.5)
            catalog['priority'][w] = 9
            catalog['code'][w] = 3

            if 'pmra' in dir(galaxy):
                #pmmem = (catalog['priority'] < 9) & ~catalog['pmra'].mask & ~catalog['pmdec'].mask
                #medianpmra = np.median(catalog['pmra'][pmmem])
                #medianpmdec = np.median(catalog['pmdec'][pmmem])
                nonmem = (catalog['code'] == 0) & (catalog['prob'] > 0) & \
                    (np.sqrt((catalog['pmra']-galaxy.pmra)**2/(catalog['pmra_error']**2 + galaxy.pmraerr**2) +
                             (catalog['pmdec']-galaxy.pmdec)**2/(catalog['pmdec_error']**2 + galaxy.pmdecerr**2)) > 3) & \
                    (catalog['pmra_error'] >= 0.0) & (catalog['pmdec_error'] >= 0.0) & ~catalog['pmra'].mask & ~catalog['pmdec'].mask
                catalog['priority'][nonmem] = 9
                        
            if showplot:    
                plt.scatter(gi0, g0, s=0.1, c=catalog['prob'], cmap='viridis')
                plt.xlabel('$(g-i)_0$')
                plt.ylabel('$g_0$')
                plt.xlim([-0.5, 2.5])
                plt.ylim([23.5, 15.5])
                plt.show()
                
                if 'pmra' in dir(galaxy):
                    w = catalog['priority'] < 9
                    plt.scatter(catalog['pmra'][w], catalog['pmdec'][w], s=0.1, c=catalog['prob'][w], cmap='viridis')
                    plt.plot(galaxy.pmra, galaxy.pmdec, 'k+')
                    plt.xlabel('PMRA')
                    plt.ylabel('PMDec')
                    plt.xlim([-5, 5])
                    plt.ylim([-5, 5])
                    plt.show()
    
        elif galaxy.name.lower() == 'ursaminor':
            predicted_y1 = -(0.65/0.6)*gi0 + 22.4
            predicted_y2 = -(1.75/0.9)*gi0 + 23.5
            predicted_y3 = -(5.75/0.3)*gi0 + 33.25
            predicted_y4 = -(5.75/0.65)*gi0 + 28.38
            
            w = (g0 > 16) & (g0 < 23) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = np.rint(11*(g0[w]-16)/(23-16)).astype(int) + 1
            
            w = (gi0 < 0.6) & ((i0 < predicted_y1) | (i0 > predicted_y2))
            catalog['priority'][w] = 13
            
            w = (gi0 >= 0.6) & ((i0 < predicted_y3) | (i0 > predicted_y4))
            catalog['priority'][w] = 13
                
            w = (g0 <= 16) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 13
            
            w = (g0 >= 23) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 13
            
            w = (catalog['cli'] > 0.5) | (catalog['clg'] > 0.5)
            catalog['priority'][w] = 14            

       
        elif galaxy.name.lower() == 'draco':
            predicted_y1 = -10*gi0 + 27.5
            predicted_y2 = -10*gi0 + 29
            
            w = (g0 > 16) & (g0 < 23) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w]=np.rint(11*(g0[w]-16)/(23-16)).astype(int) + 1
            
            w = (g0 <= 16) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 13
            
            w = (i0<predicted_y1) | (i0>predicted_y2)
            catalog['priority'][w] = 13
            
            w = (g0 >= 23) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 13
            
            w = (catalog['cli'] > 0.5) | (catalog['clg'] > 0.5)
            catalog['priority'][w] = 14
                
        elif galaxy.name.lower() == 'sculptor':
            predicted_y1 = -16*gi0 + 30.8
            predicted_y2 = -9.6*gi0 + 29.04
 
            w = (g0 > 16) & (g0 < 23) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w]=np.rint(11*(g0[w]-16)/(23-16)).astype(int) + 1
            
            w = (g0 <= 16) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 13
            
            w = (i0<predicted_y1) | (i0>predicted_y2)
            catalog['priority'][w] = 13
            
            w = (g0 >= 23) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 13
            
            w = (catalog['cli'] > 0.5) | (catalog['clg'] > 0.5)
            catalog['priority'][w] = 14
            
        elif galaxy.name.lower() == 'fornax':
            predicted_y1 = -19*gi0 + 35.15
            predicted_y2 = -9.5*gi0 + 30.875
            predicted_y3 = -3.158*gi0 + 21.68
            predicted_y4 = -1.636*gi0 + 21.045
 
            w = (g0 > 16) & (g0 < 23) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w]=np.rint(11*(g0[w]-16)/(23-16)).astype(int) + 1
            
            w = (gi0 >= 1.25) & ((i0<predicted_y3) | (i0>predicted_y4))
            catalog['priority'][w] = 13
            
            w = (g0 <= 16) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 13
            
            w = (gi0 < 1.25) & ((i0<predicted_y1) | (i0>predicted_y2))
            catalog['priority'][w] = 13
            
            w = (gi0 < 0.6) & ((i0< 20.9) | (i0>21.6))
            catalog['priority'][w] = 13
            
            w = (g0 >= 23) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 13
            
            w = (catalog['cli'] > 0.5) | (catalog['clg'] > 0.5)
            catalog['priority'][w] = 14
                
                
        elif galaxy.name.lower() == 'm31':
            catalog['prob'] = catalog['M31_RGBProb']
            top_pri = np.maximum(np.rint((i0 - 21)/(24-20) * 7).astype(int), 1) #top pri goes from 1-7 based on brightness 
            bot_pri = 12
            catalog['priority'] = bot_pri - np.rint(catalog['prob']*(bot_pri - top_pri)).astype(int)   
                        
            w = (catalog['prob'] == 0.0)
            catalog['priority'][w] = 13
            
            w = (i0 <= 20) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 13
            
            w = (i0 >= 24) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 13
            
            w = (catalog['cli'] > 0.5) | (catalog['clg'] > 0.5)
            catalog['priority'][w] = 14
            
            
        else:            
            w = (g0 > 16) & (g0 < 23) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w]=np.rint(11*(g0[w]-16)/(23-16)).astype(int) + 1
            
            w = (g0 <= 16) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 13
            
            w = (g0 >= 23) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 13
            
            w = (catalog['cli'] > 0.5) | (catalog['clg'] > 0.5)
            catalog['priority'][w] = 14

        if galaxy.name.lower() == 'fornax':
            w = (catalog['gpsf'] >= 18) & (catalog['gpsf'] <= 20) & (gi0 > 0.09) & (gi0 < 0.46) & (catalog['cli'] <= 0.5) & (catalog['clg'] <= 0.5)
            catalog['priority'][w] = 15   #ENK 9-15-2023: priority=15 is for F stars to be used as fluxstd
            
            hdul = fits.open('/home/enk/m31/dsph/for/for_moogify_member.fits.gz')
            deimos = hdul[1].data            
            c_deimos = SkyCoord(deimos['RA']*u.degree, deimos['DEC']*u.degree)
            c_pfs = SkyCoord(catalog['RA']*u.degree, catalog['Dec']*u.degree)
            idx, d2d, _ = c_deimos.match_to_catalog_sky(c_pfs)
            catalog['priority'][idx] = 0
                            
        if galaxy.name.lower() == 'ursaminor':
            hdul = fits.open('/home/enk/m31/dsph/umi/umi_moogify_member.fits.gz')
            deimos = hdul[1].data            
            c_deimos = SkyCoord(deimos['RA']*u.degree, deimos['DEC']*u.degree)
            c_pfs = SkyCoord(catalog['RA']*u.degree, catalog['Dec']*u.degree)
            idx, d2d, _ = c_deimos.match_to_catalog_sky(c_pfs)
            catalog['priority'][idx] = 0
            catalog['code'][idx] = 0
                            
        catalog['exptime'] = 1800 * np.maximum(np.minimum(np.rint(5*((i0-16) / (23.0 - 16.0)) + 1).astype(int), 6), 1)
        
        for c in catalog[(g0 < 23) & ((catalog['priority'] <= 9) & catalog['code'] == 0)]:
            s = nf.ScienceTarget(c['id'], c['RA'], c['Dec'], c['exptime'], c['priority'], 'sci')
            targets.append(s)

        if write_catalog:        
            import csv
            with open(galaxy.name.lower()+'_netflow.csv', 'w', newline='') as csvfile:
                nfwriter = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL, dialect='unix')
                nfwriter.writerow(['ob_code', 'obj_id', 'ra', 'dec', 'exptime', 'priority', 'resolution',
                                   'g_hsc', 'g_hsc_error', 'i_hsc', 'i_hsc_error', 'g_mag', 'i_mag'])
                for obj_id, c in enumerate(catalog[(g0 < 23) & (catalog['code'] == 0)]):
                    g_hsc = (c['g0']*u.ABmag).to_value(u.nJy)
                    g_hsc_error = g_hsc*c['gpsferr']*np.log(10)/2.5
                    i_hsc = (c['g0']*u.ABmag).to_value(u.nJy)
                    i_hsc_error = i_hsc*c['ipsferr']*np.log(10)/2.5
                    nfwriter.writerow([int(c['id']), obj_id+1, c['RA'], c['Dec'], c['exptime'], c['priority'], 'M',
                                      g_hsc, g_hsc_error, i_hsc, i_hsc_error, c['g0'], c['i0']])
            
        if showplot:    
            plt.scatter(gi0, g0, s=0.1, c=catalog['priority'], cmap='viridis')
            plt.xlabel('$(g-i)_0$')
            plt.ylabel('$g_0$')
            plt.xlim([-0.5, 2.5])
            plt.ylim([23.5, 15.5])
            plt.show()
    
        fluxstd = galaxy.fluxstd
        sky = galaxy.sky
        
        # class is 1 for cal class targets
        if galaxy.name.lower() == 'fornax':
            fluxstd = catalog[catalog['priority'] == 15]  #ENK 9-15-2023: pri=15 is F stars to be used as fluxstd
            fluxstd['fluxstd_id'] = fluxstd['id']
            fluxstd['ra'] = fluxstd['RA']
            fluxstd['dec'] = fluxstd['Dec']
        fluxstd['targetClass'] = 1
        fluxstd['selected'] = 0
        fluxstd['fluxstd_id'] = fluxstd['fluxstd_id'].astype(float)
        fluxstd['ra'] = fluxstd['ra'].astype(float)
        fluxstd['dec'] = fluxstd['dec'].astype(float)
        fluxstd_table = Table([fluxstd['fluxstd_id'], fluxstd['ra'], fluxstd['dec'], fluxstd['targetClass'], fluxstd['selected']],
                                 names=('id', 'RA', 'Dec', 'targetClass', 'selected'))
        removec = []
        for c in catalog.columns:
            if c.startswith('coord'):
                removec.append(c)
        catalog.remove_columns(removec)
        catalog = vstack([catalog, fluxstd_table])
        
        #class is 2 for sky class targets
        sky['targetClass'] = 2
        sky['selected'] = 0
        sky['sky_id'] = sky['sky_id'].astype(float)
        sky['ra'] = sky['ra'].astype(float)
        sky['dec'] = sky['dec'].astype(float)
        sky_table = Table([sky['sky_id'], sky['ra'], sky['dec'], sky['targetClass'], sky['selected']],
                                 names=('id', 'RA', 'Dec', 'targetClass', 'selected'))
        catalog = vstack([catalog, sky_table])
    
        #Append the fluxstd stars to the targets object
        for i in range(len(fluxstd)):
            c = nf.CalibTarget(fluxstd['fluxstd_id'][i], fluxstd['ra'][i], fluxstd['dec'][i],'cal')
            targets.append(c)
            
        #Append the sky stars to the targets object
        for i in range(len(sky)):
            c = nf.CalibTarget(sky['sky_id'][i], sky['ra'][i], sky['dec'][i],'sky') 
            targets.append(c)
    
        return(catalog, targets)
    
    
    def run_netflow(self, galaxy, catalog, targets, ra=None, dec=None, pa=0., picklefile=None, read_pickle=False):
        #Create local variables for the catalog
        #Inform the console log which galaxy is being analyzed
        if ra==None:
            ra = galaxy.ra
        if dec==None:
            dec = galaxy.dec
        
        print("The galaxy chosen is " + str(galaxy.name))
        
        bench = Bench(layout="full")
        
      
         
        # duration of one observation in seconds
        t_obs = 1800.
        
        # number of observations
        nvisit = 6
         
        telescopes = []
        forbiddenPairs = []
        for _ in range(nvisit):
            tra = ra #+ np.random.normal() * 1e-2
            tdec = dec #+ np.random.normal() * 1e-2
            posang = pa
            obs_time = "2025-05-03T10:00:00Z"
            telescope = nf.Telescope(tra, tdec, posang, obs_time)
            telescopes.append(telescope)
         
            forbiddenPairs.append([])
         
        # get focal plane positions for all targets and all visits
        fp_pos = [ t.get_fp_positions(targets) for t in telescopes ]
         
        # create the dictionary containing the costs and constraints for all classes
        # of targets
        classdict = {}
        classdict["sci_P0"] = {"nonObservationCost": 20,
                            "partialObservationCost": 1e3, "calib": False}
        classdict["sci_P1"] = {"nonObservationCost": 18,
                            "partialObservationCost": 1e3, "calib": False}
        classdict["sci_P2"] = {"nonObservationCost": 16,
                            "partialObservationCost": 1e3, "calib": False}
        classdict["sci_P3"] = {"nonObservationCost": 14,
                            "partialObservationCost": 1e3, "calib": False}
        classdict["sci_P4"] = {"nonObservationCost": 12,
                            "partialObservationCost": 1e3, "calib": False}
        classdict["sci_P5"] = {"nonObservationCost": 10,
                            "partialObservationCost": 1e3, "calib": False}
        classdict["sci_P6"] = {"nonObservationCost": 8,
                            "partialObservationCost": 1e3, "calib": False}
        classdict["sci_P7"] = {"nonObservationCost": 6,
                            "partialObservationCost": 1e3, "calib": False}
        classdict["sci_P8"] = {"nonObservationCost": 4,
                            "partialObservationCost": 1e3, "calib": False}
        classdict["sci_P9"] = {"nonObservationCost": 1,
                            "partialObservationCost": 1e3, "calib": False}
        classdict["sky"] = {"numRequired": 240,
                            "nonObservationCost": 1e4, "calib": True}
        classdict["cal"] = {"numRequired": 40,
                            "nonObservationCost": 1e4, "calib": True}
       
        # optional: slightly increase the cost for later observations,
        # to observe as early as possible
        vis_cost = [0 for i in range(nvisit)]
         
        # optional: penalize assignments where the cobra has to move far out
        def cobraMoveCost(dist):
            return 0.*dist   # Do not penalize cobra moves.
            
        gurobiOptions = dict(seed=0, presolve=1, method=4, degenmoves=0,
                            heuristics=0.5, mipfocus=3, mipgap=5.0e-04, LogToConsole=1, timelimit = 15500.0)
         
        # let's pretend that most targets have already been completely observed,
        # and that the rest has been partially observed
        alreadyObserved={}
        for t in targets:
            index = np.where(catalog['id']==t.ID)[0]  
            if index >= 0:
                alreadyObserved[t.ID] = catalog[index]['alreadyObserved']
        # for t in targets[::10]:
        #        alreadyObserved[t.ID] = 1)))

        problemfile = galaxy.name.lower()+'netflow_problem.mps'
        solutionfile = galaxy.name.lower()+'netflow_problem.sol'
         
        if read_pickle and picklefile is not None and os.path.exists(picklefile):
            with open(picklefile, 'rb') as file:
                res = pickle.load(file)
                fp_pos = pickle.load(file)
                telescopes = pickle.load(file)
        else:
            done = False
            while not done:
                # compute observation strategy
                problem = nf.buildProblem(bench, targets, fp_pos, classdict, t_obs,
                                    vis_cost, cobraMoveCost=cobraMoveCost,
                                    collision_distance=2., elbow_collisions=True,
                                    gurobi=True, gurobiOptions=gurobiOptions,
                                    alreadyObserved=None,
                                    forbiddenPairs=forbiddenPairs)
                
                print("solving the problem")
                problem.solve()
                problem._prob.write(problemfile)
                problem._prob.write(solutionfile)
                
                # extract solution
                res = [{} for _ in range(nvisit)]
                for k1, v1 in problem._vardict.items():
                    if k1.startswith("Tv_Cv_"):
                        visited = problem.value(v1) > 0
                        if visited:
                            _, _, tidx, cidx, ivis = k1.split("_")
                            res[int(ivis)][int(tidx)] = int(cidx)
                
                #catalog['id'] = catalog['id'].astype(str)
                
                done = True
    
            if picklefile is not None:
                with open(picklefile, 'wb') as file:
                    pickle.dump(res, file)
                    pickle.dump(fp_pos, file)
                    pickle.dump(telescopes, file)

        with open(galaxy.name.lower()+"_netflow_fiber_assignments.txt", "w") as f:
            print('yo')
            for i, (vis, tp, tel) in enumerate(zip(res, fp_pos, telescopes)):
                print("exposure {}:".format(i))
                print("  assigned Cobras: {}".format(len(vis)))
                tdict = defaultdict(int)
                f.write("# Exposure {}: duration {}s, RA: {}, Dec: {}, PA: {}\n".
                        format(i+1, t_obs, tel._ra, tel._dec, tel._posang))
                f.write("# Target      Fiber          X          Y         "
                        "RA        DEC\n")
                for tidx, cidx in vis.items():
                    #tdict[tgt[tidx].targetclass] += 1
                    print("{:<13n} {:6d} {:10.5f} {:10.5f} {:10.5f} {:10.5f}\n"
                            .format(int(float(targets[tidx].ID)), cidx+1, tp[tidx].real, tp[tidx].imag,
                                    targets[tidx].ra, targets[tidx].dec))
                    f.write("{:<13n} {:6d} {:10.5f} {:10.5f} {:10.5f} {:10.5f}\n"
                            .format(int(float(targets[tidx].ID)), cidx+1, tp[tidx].real, tp[tidx].imag,
                                    targets[tidx].ra, targets[tidx].dec))
                    #If the star is chosen by NetFlow, it becomes 'selected'
                    tdict[targets[tidx].targetclass] += 1
                    index = np.where(catalog['id'].astype(int)==int(float(targets[tidx].ID)))[0]
                    if len(index)>0:
                        catalog['selected'][index] = 1
                        catalog['alreadyObserved'][index] += 1
                    
                for cls, num in tdict.items():
                    print("   {}: {}".format(cls, num))
        
        return(catalog, targets)
    
    
    def count_targets(self, catalog):
        num_sky_selected = np.all((catalog['selected']==1, catalog['targetClass']==2), axis=0).sum()
        print("There are " + str(num_sky_selected) + " selected sky fibers")
        
        num_cal_selected = np.all((catalog['selected']==1, catalog['targetClass']==1), axis=0).sum()
        print("There are " + str(num_cal_selected) + " selected cal stars")
        
        count = 0
        for i in range(0,10):
            count=np.all((catalog['selected']==1, catalog['priority']==i), axis=0).sum()
            print("There are " + str(count) + " selected stars with priority " + str(i))
            count = 0
            
        count1_7 = np.all((catalog['selected']==1, catalog['targetClass']==0, catalog['priority']<=8), axis=0).sum()
        print("There are " + str(count1_7) + " selected science stars with a priority between 0 and 8")
        return()
    
    
    def plot_hexagon(self, galaxy, catalog):
        fig1, ax = plt.subplots()
        plt.xlim(galaxy.ra + np.array([-2, 2])/np.cos(np.radians(galaxy.dec)))
        plt.ylim(galaxy.dec + np.array([-2, 2]))
        ax.set_box_aspect(1)
        plt.gca().invert_xaxis()
        plt.xlabel('RA')
        plt.ylabel('Dec')
        colors = ['r', 'orange','y', 'g', 'b', 'purple','m', 'royalblue', 'c', 'pink']
        for pri in range(1,11):
            mask = catalog['priority']==pri
            plt.plot(catalog['RA'][mask], catalog['Dec'][mask], marker = '.', markersize = 0.02, color = colors[pri-1], linestyle = 'None')
        for i in range(len(galaxy.ra0)):
            thetas = np.arange(7) / 6. * 2. * np.pi
            pfs_rad = 1.44/2
            ra_hex = galaxy.ra0[i] + pfs_rad*np.sin(thetas + galaxy.pa[i]*np.pi/180) / np.cos(galaxy.dec0[i] * np.pi/180)
            dec_hex = galaxy.dec0[i] + pfs_rad*np.cos(thetas + galaxy.pa[i]*np.pi/180)
            plt.plot(ra_hex, dec_hex, 'k-', linewidth = 1)
        
        plt.show()
        
        
    def plot_netflow(self, catalog, plot_selected_only=True):
        #g0,i0 plot of science targets
        plt.xlabel(r'$(g-i)_0$')
        plt.ylabel(r'$i_0$')
        colors = ['firebrick', 'red', 'orange', 'gold', 'forestgreen', 'limegreen', 'royalblue', 'slateblue', 'blueviolet', 'gray']
        
        for pri in range(8, -1, -1):
            if plot_selected_only:
                mask = (catalog['priority']==pri) & (catalog['selected']==1)
            else:
                mask = (catalog['priority']==pri)
            i0 = np.array(catalog['i0'][mask])
            gi0 = np.array(catalog['g0'][mask] - catalog['i0'][mask])
            plt.plot(gi0, i0, color = colors[pri], marker = '.', markersize = 0.1, linestyle = 'None')
        plt.xlim([-0.5,2.5])
        plt.ylim([24,15])
        plt.show()
      
        #RA, Dec plot of science, cal, and sky targets
        plt.figure()
        plt.gca().invert_xaxis()
        plt.xlabel('RA')
        plt.ylabel('Dec')
        
        for pri in range(8, -1, -1):
            if plot_selected_only:
                mask = (catalog['priority']==pri) & (catalog['selected']==1)
            else:
                mask = (catalog['priority']==pri)
            plt.plot(catalog['RA'][mask], catalog['Dec'][mask], marker = '.', markersize = 0.1, color = colors[pri], linestyle = 'None')
        plt.show()            

        return()

GA_Netflow('ursaminor', read_pickle=False)
